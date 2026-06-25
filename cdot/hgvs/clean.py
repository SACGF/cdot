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
from dataclasses import dataclass
from enum import Enum
from typing import Optional


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
    FIXED_MULTIPLE_DOT           = "fixed_multiple_dot"
    FIXED_PREFIX_COLON           = "fixed_prefix_colon"
    STRIPPED_LEADING_JUNK        = "stripped_leading_junk"
    FIXED_GENE_WRAPPER           = "fixed_gene_wrapper"
    FIXED_SEPARATOR_TYPO         = "fixed_separator_typo"
    DROPPED_DEL_DUP_COUNT        = "dropped_del_dup_count"
    DROPPED_GENOMIC_REF_IN_PARENS = "dropped_genomic_ref_in_parens"
    ADDED_TRANSCRIPT_UNDERSCORE  = "added_transcript_underscore"
    FIXED_MULTIPLE_UNDERSCORE    = "fixed_multiple_underscore"
    UPPERCASED_BASES             = "uppercased_bases"
    RECONSTRUCTED_STRUCTURE      = "reconstructed_structure"
    ADDED_MISSING_KIND           = "added_missing_kind"
    # Gene/transcript fixes (all WARNING)
    SWAPPED_GENE_TRANSCRIPT      = "swapped_gene_transcript"
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
    # Gene → transcript resolution (WARNING on success, ERROR on failure)
    RESOLVED_GENE_TO_TRANSCRIPT  = "resolved_gene_to_transcript"
    NO_TRANSCRIPT_FOR_GENE       = "no_transcript_for_gene"
    # Validation errors (ERROR)
    NO_COLON                     = "no_colon"
    MISSING_REFERENCE_SEQUENCE   = "missing_reference_sequence"
    INS_WITH_INTEGER_LENGTH      = "ins_with_integer_length"
    INS_MISSING_SEQUENCE         = "ins_missing_sequence"


class HGVSCleanOp(Enum):
    """A single cleaning operation in the clean_hgvs() pipeline.

    Unlike HGVSFixCode (which describes *outcomes* and is not 1:1 with steps),
    each member here names exactly one step, so it is a stable key for selecting
    a subset of cleaning operations.  Members are listed in pipeline order, but
    selection only *filters* the pipeline — it never reorders it (several steps
    depend on earlier ones running first).
    """
    STRIP_PROTEIN_SUFFIX     = "strip_protein_suffix"
    STRIP_LEADING_JUNK       = "strip_leading_junk"
    STRIP_WHITESPACE         = "strip_whitespace"
    STRIP_SURROUNDING_PUNCTUATION = "strip_surrounding_punctuation"
    FIX_MULTIPLE_COLON       = "fix_multiple_colon"
    FIX_GENE_WRAPPER         = "fix_gene_wrapper"
    FIX_MULTIPLE_DOT         = "fix_multiple_dot"
    FIX_MULTIPLE_UNDERSCORE  = "fix_multiple_underscore"
    FIX_PREFIX_COLON         = "fix_prefix_colon"
    FIX_MULTIPLE_KIND        = "fix_multiple_kind"
    FIX_SEPARATOR_TYPO       = "fix_separator_typo"
    ADD_N_PREFIX             = "add_n_prefix"
    LOWERCASE_MUTATION_TYPE  = "lowercase_mutation_type"
    DROP_DEL_DUP_COUNT       = "drop_del_dup_count"
    STRIP_UNBALANCED_BRACKETS = "strip_unbalanced_brackets"
    ADD_TRANSCRIPT_UNDERSCORE = "add_transcript_underscore"
    RECONSTRUCT_STRUCTURE    = "reconstruct_structure"
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


def _looks_like_transcript(s: str) -> bool:
    """True if s starts with a known mRNA/RNA RefSeq prefix, or is Ensembl/LRG.

    Case-sensitive — callers must uppercase s themselves if they want a
    case-insensitive check (matching VG behaviour).
    """
    if s.startswith('ENST') or s.startswith('LRG_'):
        return True
    # Keys in _REFSEQ_PREFIX_TYPE are uppercase 3-char prefixes e.g. 'NM_'
    return _REFSEQ_PREFIX_TYPE.get(s[:3]) in ('mRNA', 'RNA')


def _looks_like_hgvs_prefix(s: str) -> bool:
    """True if s starts with any known RefSeq prefix, or is Ensembl/LRG.

    Case-sensitive — callers must uppercase s themselves if needed.
    """
    if s.startswith('ENST') or s.startswith('LRG_'):
        return True
    return s[:3] in _REFSEQ_PREFIX_TYPE


# ---------------------------------------------------------------------------
# Compiled regex patterns
# ---------------------------------------------------------------------------

_P_HGVS_REMOVAL = re.compile(r"^(?P<main>.*?) (p\.|[(]p\.).*$")

# Two-letter RefSeq accession prefixes (NM, NR, XM, XR, NC, NG, ...) taken from
# the type table, used to repair an accession missing its underscore. Requires
# 6+ digits (every RefSeq accession has at least six) so a real gene symbol that
# starts with these letters plus a stray digit (eg NR3C1, NS1, ZP3) is never
# mistaken for an accession and mangled.
_TWO_LETTER_PREFIXES = sorted({p[:2] for p in _REFSEQ_PREFIX_TYPE})
_TRANSCRIPT_NO_UNDERSCORE = re.compile(
    rf"^({'|'.join(_TWO_LETTER_PREFIXES)})(\d{{6,}})", re.IGNORECASE)

# Second letters for which prepending "N" yields a real genomic/transcript prefix
# (NM, NC, NR, NG, NT, NW, NS, NZ, NP) — used to restore a dropped leading "N".
_N_PREFIX_SECOND_LETTERS = frozenset(
    p[1] for p in _REFSEQ_PREFIX_TYPE if p.startswith("N"))

# Lenient pattern — matches HGVS strings whether or not they have a transcript,
# gene symbol, version, colon, or dot.  Used to reconstruct the canonical form.
_HGVS_CLEAN_PATTERN = re.compile(
    # transcript_value allows an optional LRG transcript/protein suffix (eg the
    # "t1" in "LRG_199t1") so the suffix isn't mis-parsed as the gene symbol.
    r"^((?P<transcript_prefix>NR|NM|NC|ENST|LRG|XR|XM)(?P<transcript_value>[0-9_]*(?:[tp][0-9]+)?)"
    r"(\.(?P<transcript_version>[0-9]+))?)?"
    r"([(]?(?P<gene_symbol>[^):\.]{2,8})[)]?)?"
    # A stray colon between the kind letter and the coordinate (eg "...2c:1188")
    # is accepted as the kind separator, so the same glued-kind input reconstructs
    # whether or not that colon is present (it is just misplaced, belonging before
    # the kind letter rather than after it).
    r":?(?P<letter>c|g|n|p)[.:]?(?P<nomen>[^:\.]*)$",
    re.IGNORECASE,
)

# Missing kind letter after the colon, classified by reference-sequence type so
# the right kind is inserted: mRNA/Ensembl transcript → "c.", non-coding RNA →
# "n.", genomic → "g.". The coordinate may start with a digit, "-" (5' UTR) or
# "*" (3' UTR).
_COORD_START = r"[\d*\-]"
_HGVS_MRNA_NO_KIND    = re.compile(rf"^(NM_|XM_|ENST)\d+.*:{_COORD_START}", re.IGNORECASE)
_HGVS_RNA_NO_KIND     = re.compile(rf"^(NR_|XR_)\d+.*:{_COORD_START}", re.IGNORECASE)
_HGVS_GENOMIC_NO_KIND = re.compile(
    rf"^(NC_|NG_|NW_|NT_|AC_|NS_|NZ_)\d+.*:{_COORD_START}", re.IGNORECASE)

_REF_ALT_NUC     = re.compile(r"(?P<ref>[gatc]+)>(?P<alt>[gatc=]+)$", re.IGNORECASE)
_DEL_INS_DUP_NUC = re.compile(
    r"(del(?P<del>[gatc]+))?(?P<op>ins|dup|del)(?P<ins>[gatc]*)$",
    re.IGNORECASE,
)

_GENE_IN_PARENS = re.compile(r"^(.+)\((.*)\)$")

_INS_INTEGER = re.compile(r".*ins\d+$", re.IGNORECASE)
_INS_EMPTY   = re.compile(r".*ins$",    re.IGNORECASE)

# Misplaced colon between an mRNA/RNA prefix and its underscore, e.g. "NM:_000059"
_PREFIX_COLON_UNDERSCORE = re.compile(r"^(N[MR]|X[MR]):_", re.IGNORECASE)

# Transcript accession prefixes (used to locate the HGVS core within a string)
_TX_PREFIX = r"(?:NM_|NR_|XM_|XR_|NC_|NG_|ENST|LRG_)"

# Leading junk before a transcript accession, e.g. "#", ":", a literal "\t",
# or a genome-build prefix like "GRCh38.p2 " / "GRCh37.p12#"
_LEADING_JUNK = re.compile(r"^(?:GRCh\d+(?:\.p\d+)?[#\s]*|\\t|[#:])+", re.IGNORECASE)

# Gene symbol wedged between extra colons, e.g. "tx:(GENE):c." / "tx:(GENE)c." /
# "tx:GENE:c." — should be "tx(GENE):c."
_GENE_WRAPPER = re.compile(
    rf"^({_TX_PREFIX}[\w.]+):\(?([A-Za-z][A-Za-z0-9-]+)\)?:?([cgnmp][.,])", re.IGNORECASE)
# Stray colon inside the gene parens, e.g. "(RAD51C:)" -> "(RAD51C)"
_GENE_PARENS_COLON = re.compile(r"\(([A-Za-z][A-Za-z0-9-]*):\)")
# Parens wrapping the whole accession, e.g. "(NM_000059.4):c." -> "NM_000059.4:c."
_PARENS_WRAPPED_ACCESSION = re.compile(rf"^\(({_TX_PREFIX}[\w.]+)\)(?=:)", re.IGNORECASE)
# Missing colon between a (gene) wrapper and the kind, e.g.
# "BRCA1(NM_000059.4)c.123A>G" -> "BRCA1(NM_000059.4):c.123A>G" (the
# gene/transcript swap is then done by _fix_gene_transcript).
_MISSING_COLON_BEFORE_KIND = re.compile(r"(\([^)]+\))([cgnmp]\.)", re.IGNORECASE)

# Comma used in place of the c./g. dot, e.g. ":c,1811" -> ":c.1811"
_KIND_COMMA = re.compile(r":([cgnmp]),", re.IGNORECASE)
# Colon used in place of the kind dot, e.g. ":c:1811" -> ":c.1811". The
# lookahead keeps this from touching a normal ":c." (the second char is a dot,
# not a coordinate).
_KIND_COLON = re.compile(r":([cgnmp]):(?=[\d*?(\-])", re.IGNORECASE)
# Period used in place of the substitution '>', e.g. "1030C.T" -> "1030C>T"
_SUB_PERIOD = re.compile(r"([ACGT])\.([ACGT])$", re.IGNORECASE)

# del/dup with an explicit range followed by a redundant base count, e.g.
# "c.1315_1337dup23" -> "c.1315_1337dup" (the range already gives the length;
# normalise to plain del/dup). The range guard avoids changing the meaning of a
# single-position "c.123del5" (delete 5 bases). delins/ins are excluded — they
# genuinely need the inserted sequence, so they stay flagged as errors.
_DEL_DUP_RANGE_COUNT = re.compile(r"(\d+_\d[\d+\-_*]*(?:del|dup))\d+$", re.IGNORECASE)

# A genomic accession (NC_/NG_/NW_) wedged into the gene-symbol parenthetical
# slot of a transcript reference, e.g. "NM_000059.4(NC_000013.11):c.68del". A
# genomic accession is not a gene symbol, so the parser rejects it; dropping the
# parenthetical (the transcript alone is a complete reference) makes it parse.
# Only fires when the part *before* the parens is a transcript and the part
# *inside* is genomic — this deliberately does NOT touch the valid HGVS form
# "NC_000013.11(NM_000059.4):c.68del" (genomic reference, transcript selector),
# where the parenthetical holds the transcript, not the genomic accession.
_GENOMIC_REF_IN_GENE_PARENS = re.compile(
    r"^((?:NM_|NR_|XM_|XR_|ENST|LRG_)[\w.]+)\((?:NC_|NG_|NW_)\d[\w.]*\)(?=:)",
    re.IGNORECASE,
)


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


def _reconstruct(m: re.Match) -> str:
    """Rebuild a canonical HGVS string from a _HGVS_CLEAN_PATTERN match."""
    kind = m.group("letter").lower()
    nomen = m.group("nomen")
    c_nomen = f":{kind}.{nomen}"

    gene_symbol = m.group("gene_symbol")
    transcript_prefix = m.group("transcript_prefix")

    if transcript_prefix:
        transcript_value = m.group("transcript_value")
        text = f"{transcript_prefix.upper()}{transcript_value}"
        if transcript_version := m.group("transcript_version"):
            text += f".{transcript_version}"
        if gene_symbol:
            text += f"({gene_symbol})"
        return text + c_nomen
    elif gene_symbol:
        return f"{gene_symbol}{c_nomen}"
    return None  # can't reconstruct without transcript or gene


def _fix_gene_transcript(hgvs_string: str) -> tuple[str, list[HGVSFix]]:
    """
    Fix the common clinical mistake of swapping gene symbol and transcript,
    and fix lowercase transcript identifiers.

    Ported from HGVSMatcher.fix_gene_transcript() in vg_code/hgvs/hgvs_matcher.py.
    """
    fixes = []

    try:
        prefix, allele = hgvs_string.split(":", 1)
    except ValueError:
        return hgvs_string, []

    if m := _GENE_IN_PARENS.match(prefix):
        transcript_accession, gene_symbol = m.groups()
    else:
        transcript_accession, gene_symbol = prefix, None

    transcript_ok = _looks_like_transcript(transcript_accession)

    if not transcript_ok and gene_symbol:
        uc_gene = gene_symbol.upper()
        if _looks_like_transcript(uc_gene):
            # Gene and transcript are swapped
            transcript_accession, gene_symbol = gene_symbol, transcript_accession
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=HGVSFixCode.SWAPPED_GENE_TRANSCRIPT,
                message="Swapped gene/transcript",
                original=prefix,
                fixed=f"{transcript_accession}({gene_symbol})" if gene_symbol else transcript_accession,
            ))
            transcript_ok = _looks_like_transcript(transcript_accession)
        elif _looks_like_hgvs_prefix(uc_gene):
            gene_symbol = uc_gene
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=HGVSFixCode.UPPERCASED_HGVS_PREFIX,
                message="Upper cased HGVS prefix",
                original=prefix,
                fixed=f"{transcript_accession}({gene_symbol})",
            ))

    if not transcript_ok and transcript_accession:
        uc = transcript_accession.upper()
        if _looks_like_transcript(uc):
            transcript_accession = uc
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=HGVSFixCode.UPPERCASED_TRANSCRIPT,
                message="Upper cased transcript",
                original=prefix,
                fixed=transcript_accession,
            ))

    new_prefix = f"{transcript_accession}({gene_symbol})" if gene_symbol else transcript_accession
    return f"{new_prefix}:{allele}", fixes


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
    for ch in (":", "c", "g", "n", "m", "p", "r", "."):
        if s.startswith(ch):
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.ERROR,
                code=HGVSFixCode.MISSING_REFERENCE_SEQUENCE,
                message=f"Missing reference sequence (string starts with '{ch}')",
            ))
            break
    return fixes


# ---------------------------------------------------------------------------
# Cleaning operations
#
# Each op has the signature (s: str) -> tuple[str, list[HGVSFix]]: it returns
# the (possibly unchanged) string and a list of HGVSFix describing what it did
# (empty if it made no change).  Ops are pure and order-dependent — _PIPELINE
# fixes the canonical order; selection only filters which ops run.
# ---------------------------------------------------------------------------

def _op_strip_protein_suffix(s: str) -> tuple[str, list[HGVSFix]]:
    # e.g. "NM_000059.4:c.316+5G>A p.Arg106*" → "NM_000059.4:c.316+5G>A"
    if m := _P_HGVS_REMOVAL.match(s):
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
    # "NM_000059.4:c.123A>G;". No valid HGVS string starts/ends with these.
    s2 = s.strip().strip("'\"`").rstrip("/;,").strip()
    if s2 != s:
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
    for kind in "cgmnpr":
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
    return s, fixes


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


def _op_lowercase_mutation_type(s: str) -> tuple[str, list[HGVSFix]]:
    # Lowercase mutation types (DUP→dup, DEL→del, INS→ins, INV→inv)
    # Also handles DELINS since it replaces substrings.
    # Only operate on the allele (the part after the first ':'): gene symbols such
    # as INSR, INVS, INSL3 or DELEC1 contain a mutation-type substring and must
    # not be corrupted (e.g. "NM_000208.4(INSR):c.215A>G" must not become "insR").
    if ":" not in s:
        return s, []
    prefix, allele = s.split(":", 1)
    fixes = []
    for mt in ["INS", "DEL", "DUP", "INV"]:
        if mt in allele:
            allele2 = allele.replace(mt, mt.lower())
            new_s = f"{prefix}:{allele2}"
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=HGVSFixCode.LOWERCASED_MUTATION_TYPE,
                message=f"Lowercased mutation type '{mt}'",
                original=s, fixed=new_s,
            ))
            allele = allele2
            s = new_s
    return s, fixes


def _op_strip_unbalanced_brackets(s: str) -> tuple[str, list[HGVSFix]]:
    # Only strip when brackets are genuinely unbalanced (a stray "(" or ")" from
    # copy/paste). Do NOT strip when they balance: multiple balanced pairs are
    # valid HGVS - a gene symbol in parens plus uncertain-range notation, e.g.
    # "NM_004006.2(DMD):c.(4071+1_4072-1)_(5154+1_5155-1)del".
    open_b  = s.count("(")
    close_b = s.count(")")
    if open_b != close_b:
        s2 = s.replace("(", "").replace(")", "")
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.STRIPPED_UNBALANCED_BRACKETS,
            message="Stripped unbalanced brackets",
            original=s, fixed=s2,
        )]
    return s, []


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


def _op_reconstruct_structure(s: str) -> tuple[str, list[HGVSFix]]:
    # Structural reconstruction: uppercase prefix, lowercase kind letter,
    # insert ':' and '.' where missing, place gene symbol correctly
    if m := _HGVS_CLEAN_PATTERN.match(s):
        # A 2-letter RefSeq prefix (NM, NR, NC, XM, ...) is always followed by an
        # underscore, and _op_add_transcript_underscore has already inserted any
        # missing one by now. If the matched "transcript" still has none, this is
        # really a gene symbol that merely starts with those letters plus a digit
        # (eg NR3C1, NC2) — do not reconstruct it into "NR3(C1)".
        prefix = m.group("transcript_prefix")
        value = m.group("transcript_value") or ""
        if prefix and len(prefix) == 2 and "_" not in value:
            return s, []
        reconstructed = _reconstruct(m)
        if reconstructed and reconstructed != s:
            transcript_prefix = m.group("transcript_prefix")
            if transcript_prefix and transcript_prefix != transcript_prefix.upper():
                code = HGVSFixCode.UPPERCASED_TRANSCRIPT
                message = "Uppercased transcript prefix"
            else:
                code = HGVSFixCode.RECONSTRUCTED_STRUCTURE
                message = "Reconstructed HGVS with correct casing and structure"
            return reconstructed, [HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=code, message=message,
                original=s, fixed=reconstructed,
            )]
    return s, []


def _op_uppercase_bases(s: str) -> tuple[str, list[HGVSFix]]:
    fixes = []
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


def _op_add_missing_kind(s: str) -> tuple[str, list[HGVSFix]]:
    # Insert a missing kind letter + dot after the colon, choosing the kind from
    # the reference-sequence type:
    #   mRNA / Ensembl transcript (NM_, XM_, ENST) → "c."  eg "NM_001754.5:557T>A"
    #   non-coding RNA           (NR_, XR_)        → "n."  eg "NR_002196.2:601A>G"
    #   genomic                  (NC_, NG_, ...)   → "g."  eg "NC_000007.13:1G>T"
    for pattern, kind in (
        (_HGVS_MRNA_NO_KIND, "c"),
        (_HGVS_RNA_NO_KIND, "n"),
        (_HGVS_GENOMIC_NO_KIND, "g"),
    ):
        if pattern.match(s):
            s2 = s.replace(":", f":{kind}.", 1)
            return s2, [HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=HGVSFixCode.ADDED_MISSING_KIND,
                message=f"Added missing '{kind}.' kind prefix",
                original=s, fixed=s2,
            )]
    return s, []


def _op_strip_leading_junk(s: str) -> tuple[str, list[HGVSFix]]:
    # Strip junk before the transcript accession, e.g. "#NM_...", ":NM_...",
    # "GRCh38.p2 NM_...". Only fires when a transcript prefix actually follows,
    # so a leading gene symbol (e.g. "BRCA1(NM_...)") is never touched.
    s2 = _LEADING_JUNK.sub("", s)
    if s2 != s and re.match(_TX_PREFIX, s2, re.IGNORECASE):
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
    #   "(NM_000059.4):c."                           -> "NM_000059.4:c."
    #   "BRCA1(NM_000059.4)c."                       -> "BRCA1(NM_000059.4):c."
    s2 = _GENE_PARENS_COLON.sub(r"(\1)", s)
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


def _op_fix_separator_typo(s: str) -> tuple[str, list[HGVSFix]]:
    # Comma ("c,1811") or colon ("c:1811") in place of the kind dot, and period
    # in place of the substitution '>' ("1030C.T" -> "1030C>T").
    s2 = _KIND_COMMA.sub(r":\1.", s)
    s2 = _KIND_COLON.sub(r":\1.", s2)
    s2 = _SUB_PERIOD.sub(r"\1>\2", s2)
    if s2 != s:
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.FIXED_SEPARATOR_TYPO,
            message="Fixed comma/colon/period separator typo", original=s, fixed=s2,
        )]
    return s, []


def _op_drop_del_dup_count(s: str) -> tuple[str, list[HGVSFix]]:
    # Drop a redundant base count after a ranged del/dup, normalising to plain
    # del/dup, e.g. "c.1315_1337dup23" -> "c.1315_1337dup".
    head, sep, allele = s.rpartition(":")
    target = allele if sep else s
    if _DEL_DUP_RANGE_COUNT.search(target):
        new_allele = _DEL_DUP_RANGE_COUNT.sub(r"\1", target)
        s2 = f"{head}:{new_allele}" if sep else new_allele
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.DROPPED_DEL_DUP_COUNT,
            message="Dropped redundant base count after ranged del/dup",
            original=s, fixed=s2,
        )]
    return s, []


def _op_drop_genomic_ref_in_parens(s: str) -> tuple[str, list[HGVSFix]]:
    # Drop a genomic accession sitting in the gene-symbol parenthetical slot,
    # e.g. "NM_000059.4(NC_000013.11):c.68del" -> "NM_000059.4:c.68del". The
    # transcript on its own is a complete reference, so the genomic accession is
    # redundant and (being neither a gene symbol nor a valid selector here) only
    # blocks parsing.
    s2 = _GENOMIC_REF_IN_GENE_PARENS.sub(r"\1", s)
    if s2 != s:
        return s2, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.DROPPED_GENOMIC_REF_IN_PARENS,
            message="Dropped genomic accession from the gene-symbol parentheses",
            original=s, fixed=s2,
        )]
    return s, []


# Canonical pipeline: (op, function), in the order they must run.
# clean_hgvs() filters this by the caller's `ops` set but never reorders it.
_PIPELINE: list[tuple[HGVSCleanOp, "callable"]] = [
    (HGVSCleanOp.STRIP_PROTEIN_SUFFIX,      _op_strip_protein_suffix),
    (HGVSCleanOp.STRIP_LEADING_JUNK,        _op_strip_leading_junk),
    (HGVSCleanOp.STRIP_WHITESPACE,          _op_strip_whitespace),
    (HGVSCleanOp.STRIP_SURROUNDING_PUNCTUATION, _op_strip_surrounding_punctuation),
    (HGVSCleanOp.FIX_MULTIPLE_COLON,        _op_fix_multiple_colon),
    (HGVSCleanOp.FIX_GENE_WRAPPER,          _op_fix_gene_wrapper),
    (HGVSCleanOp.FIX_MULTIPLE_DOT,          _op_fix_multiple_dot),
    (HGVSCleanOp.FIX_MULTIPLE_UNDERSCORE,   _op_fix_multiple_underscore),
    (HGVSCleanOp.FIX_PREFIX_COLON,          _op_fix_prefix_colon),
    (HGVSCleanOp.FIX_MULTIPLE_KIND,         _op_fix_multiple_kind),
    (HGVSCleanOp.FIX_SEPARATOR_TYPO,        _op_fix_separator_typo),
    (HGVSCleanOp.ADD_N_PREFIX,              _op_add_n_prefix),
    (HGVSCleanOp.LOWERCASE_MUTATION_TYPE,   _op_lowercase_mutation_type),
    (HGVSCleanOp.DROP_DEL_DUP_COUNT,        _op_drop_del_dup_count),
    (HGVSCleanOp.STRIP_UNBALANCED_BRACKETS, _op_strip_unbalanced_brackets),
    (HGVSCleanOp.ADD_TRANSCRIPT_UNDERSCORE, _op_add_transcript_underscore),
    (HGVSCleanOp.RECONSTRUCT_STRUCTURE,     _op_reconstruct_structure),
    (HGVSCleanOp.UPPERCASE_BASES,           _op_uppercase_bases),
    (HGVSCleanOp.ADD_MISSING_KIND,          _op_add_missing_kind),
    (HGVSCleanOp.FIX_GENE_TRANSCRIPT,       _fix_gene_transcript),
    # Runs after the gene/transcript swap: an "NC_…(NM_…)" genomic-reference form
    # is normalised to "NM_…(NC_…)" by that step, so dropping the genomic
    # parenthetical here collapses both the broken and the swapped form to the
    # same parseable transcript reference.
    (HGVSCleanOp.DROP_GENOMIC_REF_IN_PARENS, _op_drop_genomic_ref_in_parens),
]

# All cleaning ops — the default for clean_hgvs(). Use set algebra for a
# blocklist, e.g. ALL_CLEAN_OPS - {HGVSCleanOp.STRIP_PROTEIN_SUFFIX}.
ALL_CLEAN_OPS: frozenset = frozenset(HGVSCleanOp)


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
            run in the canonical pipeline order regardless of set order.
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

    fixes: list[HGVSFix] = []
    s = hgvs_string

    for op, fn in _PIPELINE:
        if op in ops:
            s, op_fixes = fn(s)
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
    layer — ``resolve_transcript_version`` calls the data provider's
    ``is_version_substitution_safe`` and upgrades this plain USED_ADJACENT_VERSION
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
