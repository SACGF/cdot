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
    FIXED_DOUBLE_COLON           = "fixed_double_colon"
    FIXED_DOUBLE_KIND            = "fixed_double_kind"
    ADDED_N_PREFIX               = "added_n_prefix"
    LOWERCASED_MUTATION_TYPE     = "lowercased_mutation_type"
    STRIPPED_UNBALANCED_BRACKETS = "stripped_unbalanced_brackets"
    ADDED_TRANSCRIPT_UNDERSCORE  = "added_transcript_underscore"
    UPPERCASED_BASES             = "uppercased_bases"
    ADDED_MISSING_KIND           = "added_missing_kind"
    # Gene/transcript fixes (all WARNING)
    SWAPPED_GENE_TRANSCRIPT      = "swapped_gene_transcript"
    UPPERCASED_TRANSCRIPT        = "uppercased_transcript"
    UPPERCASED_HGVS_PREFIX       = "uppercased_hgvs_prefix"
    # Transcript version fallback (WARNING)
    USED_ADJACENT_VERSION        = "used_adjacent_version"
    # Gene → transcript resolution (WARNING on success, ERROR on failure)
    RESOLVED_GENE_TO_TRANSCRIPT  = "resolved_gene_to_transcript"
    NO_TRANSCRIPT_FOR_GENE       = "no_transcript_for_gene"
    # Validation errors (ERROR)
    NO_COLON                     = "no_colon"
    MISSING_REFERENCE_SEQUENCE   = "missing_reference_sequence"
    INS_WITH_INTEGER_LENGTH      = "ins_with_integer_length"
    INS_MISSING_SEQUENCE         = "ins_missing_sequence"


@dataclass(frozen=True)
class HGVSFix:
    severity: HGVSFixSeverity
    code: HGVSFixCode
    message: str
    original: Optional[str] = None  # before-value, if applicable
    fixed: Optional[str] = None     # after-value, if applicable


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

_TRANSCRIPT_NO_UNDERSCORE = re.compile(r"^(NM|NC)(\d+)", re.IGNORECASE)

# Lenient pattern — matches HGVS strings whether or not they have a transcript,
# gene symbol, version, colon, or dot.  Used to reconstruct the canonical form.
_HGVS_CLEAN_PATTERN = re.compile(
    r"^((?P<transcript_prefix>NR|NM|NC|ENST|LRG|XR)(?P<transcript_value>[0-9_]*)"
    r"(\.(?P<transcript_version>[0-9]+))?)?"
    r"([(]?(?P<gene_symbol>[^):\.]{2,8})[)]?)?"
    r":?(?P<letter>c|g|n|p)[.]?(?P<nomen>[^:\.]*)$",
    re.IGNORECASE,
)

_HGVS_TRANSCRIPT_NO_CDOT = re.compile(r"^(NM_|ENST)\d+.*:\d+")
_HGVS_CONTIG_NO_GDOT      = re.compile(r"^NC_\d+.*:\d+")

_REF_ALT_NUC     = re.compile(r"(?P<ref>[gatc]+)>(?P<alt>[gatc=]+)$", re.IGNORECASE)
_DEL_INS_DUP_NUC = re.compile(
    r"(del(?P<del>[gatc]+))?(?P<op>ins|dup|del)(?P<ins>[gatc]*)$",
    re.IGNORECASE,
)

_GENE_IN_PARENS = re.compile(r"^(.+)\((.*)\)$")

_INS_INTEGER = re.compile(r".*ins\d+$", re.IGNORECASE)
_INS_EMPTY   = re.compile(r".*ins$",    re.IGNORECASE)


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
    for ch in (":", "c", "g", "."):
        if s.startswith(ch):
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.ERROR,
                code=HGVSFixCode.MISSING_REFERENCE_SEQUENCE,
                message=f"Missing reference sequence (string starts with '{ch}')",
            ))
            break
    return fixes


# ---------------------------------------------------------------------------
# Public: clean_hgvs
# ---------------------------------------------------------------------------

def clean_hgvs(
    hgvs_string: str,
    raise_on_errors: bool = False,
) -> tuple[str, list[HGVSFix]]:
    """
    Attempt to fix common formatting problems in an HGVS string.

    Returns:
        (cleaned_string, fixes)
        where fixes is a list of HGVSFix describing every change made and
        any validation errors found.

    If raise_on_errors=True, raises HGVSInputError on the first ERROR-level
    fix instead of returning it in the list.

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
    fixes: list[HGVSFix] = []
    s = hgvs_string

    # 1. Strip trailing protein (p.) annotation separated by a space
    #    e.g. "NM_000059.4:c.316+5G>A p.Arg106*" → "NM_000059.4:c.316+5G>A"
    if m := _P_HGVS_REMOVAL.match(s):
        new_s = m.group("main")
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.STRIPPED_PROTEIN_SUFFIX,
            message="Removed trailing protein (p.) annotation",
            original=s, fixed=new_s,
        ))
        s = new_s

    # 2. Strip non-printable characters and all whitespace
    s2 = "".join(c for c in s if c.isprintable()).replace(" ", "")
    if s2 != s:
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.STRIPPED_WHITESPACE,
            message="Removed whitespace",
            original=s, fixed=s2,
        ))
        s = s2

    # 3. Fix "::" → ":"
    if "::" in s:
        s2 = s.replace("::", ":")
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.FIXED_DOUBLE_COLON,
            message="Fixed double colon '::'", original=s, fixed=s2,
        ))
        s = s2

    # 3b. Fix "__" → "_"
    if "__" in s:
        s2 = s.replace("__", "_")
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.ADDED_TRANSCRIPT_UNDERSCORE,
            message="Fixed double underscore '__'", original=s, fixed=s2,
        ))
        s = s2

    # 4. Fix double-kind notation (e.g. "c.c." → "c.")
    for kind in "cgmnpr":
        pat = f"{kind}.{kind}."
        if pat in s.lower():
            s2 = re.sub(re.escape(pat), f"{kind}.", s, flags=re.IGNORECASE)
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=HGVSFixCode.FIXED_DOUBLE_KIND,
                message=f"Fixed double kind notation '{pat}'",
                original=s, fixed=s2,
            ))
            s = s2

    # 5. Add missing "N" prefix (e.g. "M_001234" → "NM_001234")
    if s[:2].upper() in ("M_", "C_", "R_"):
        s2 = "N" + s
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.ADDED_N_PREFIX,
            message=f"Added missing 'N' prefix",
            original=s, fixed=s2,
        ))
        s = s2

    # 6. Lowercase mutation types (DUP→dup, DEL→del, INS→ins, INV→inv)
    #    Also handles DELINS since it replaces substrings
    for mt in ["INS", "DEL", "DUP", "INV"]:
        if mt in s:
            s2 = s.replace(mt, mt.lower())
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=HGVSFixCode.LOWERCASED_MUTATION_TYPE,
                message=f"Lowercased mutation type '{mt}'",
                original=s, fixed=s2,
            ))
            s = s2

    # 7. Strip unbalanced or multiple brackets
    open_b  = s.count("(")
    close_b = s.count(")")
    if open_b != close_b or open_b > 1 or close_b > 1:
        s2 = s.replace("(", "").replace(")", "")
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.STRIPPED_UNBALANCED_BRACKETS,
            message="Stripped unbalanced or multiple brackets",
            original=s, fixed=s2,
        ))
        s = s2

    # 8. Add missing underscore to NM/NC prefix (e.g. "NM12345" → "NM_12345")
    s2 = _TRANSCRIPT_NO_UNDERSCORE.sub(r"\1_\2", s)
    if s2 != s:
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.ADDED_TRANSCRIPT_UNDERSCORE,
            message="Added missing underscore to transcript prefix",
            original=s, fixed=s2,
        ))
        s = s2

    # 9. Structural reconstruction: uppercase prefix, lowercase kind letter,
    #    insert ':' and '.' where missing, place gene symbol correctly
    if m := _HGVS_CLEAN_PATTERN.match(s):
        reconstructed = _reconstruct(m)
        if reconstructed and reconstructed != s:
            transcript_prefix = m.group("transcript_prefix")
            if transcript_prefix and transcript_prefix != transcript_prefix.upper():
                code = HGVSFixCode.UPPERCASED_TRANSCRIPT
                message = "Uppercased transcript prefix"
            else:
                code = HGVSFixCode.UPPERCASED_BASES
                message = "Reconstructed HGVS with correct casing and structure"
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=code, message=message,
                original=s, fixed=reconstructed,
            ))
            s = reconstructed

    # 10. Uppercase bases in ref>alt substitutions (e.g. "c.100a>g" → "c.100A>G")
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

    # 11. Insert missing "c." or "g." after the colon
    #     e.g. "NM_001754.5:557T>A" → "NM_001754.5:c.557T>A"
    if _HGVS_TRANSCRIPT_NO_CDOT.match(s):
        s2 = s.replace(":", ":c.", 1)
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.ADDED_MISSING_KIND,
            message="Added missing 'c.' kind prefix",
            original=s, fixed=s2,
        ))
        s = s2
    elif _HGVS_CONTIG_NO_GDOT.match(s):
        s2 = s.replace(":", ":g.", 1)
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.ADDED_MISSING_KIND,
            message="Added missing 'g.' kind prefix",
            original=s, fixed=s2,
        ))
        s = s2

    # 12. Gene/transcript swap and case normalisation
    s, gene_fixes = _fix_gene_transcript(s)
    fixes.extend(gene_fixes)

    # 13. Validation — run after cleaning so both the fix and any remaining
    #     problem are captured in one call
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


def _rank_versions(
    requested: int,
    available: list[int],
    strategy: VersionStrategy,
) -> list[int]:
    """
    Return available versions sorted by preference given the strategy.

    Ported from HGVSMatcher._get_sort_key_transcript_version_and_methods()
    in vg_code/hgvs/hgvs_matcher.py (minus the method/ClinGen sorting,
    which is VG-specific).
    """
    if strategy == VersionStrategy.LATEST:
        return sorted(available, reverse=True)

    def sort_key(v: int) -> tuple:
        distance = abs(requested - v)
        prefer_later = v < requested  # False=later, True=earlier → sort asc = later first
        if strategy == VersionStrategy.CLOSEST:
            return (distance, prefer_later)
        else:  # UP_THEN_DOWN
            return (prefer_later, distance)

    return sorted(available, key=sort_key)


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

    ranked = _rank_versions(requested_version, available_versions, strategy)
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
