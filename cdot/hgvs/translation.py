"""
c. to p. conversion that uses cdot's translation data (data schema >= 0.2.35).

biocommons HGVS translates the CDS of the transcript sequence with one genetic code, and
doesn't know about RefSeq translation exceptions. fix_c_to_p() wraps a biocommons VariantMapper
or AssemblyMapper c_to_p() call, and uses cdot's transcript data to:

* pick the genetic code: the vertebrate mitochondrial code for ``transl_table`` 2, and the
  selenocysteine code for selenoproteins (biocommons only detects these when it can fetch the
  protein sequence, which eg FastaSeqFetcher can't provide, and otherwise reads Sec as a stop)
* report a new UGA in a selenoprotein as a stop, as it is only read as Sec at the known codons
* refuse transcripts with a ribosomal frameshift, which biocommons can't translate
* check the codons RefSeq gives a translation exception for against the transcript sequence
  used. A genome derived sequence (eg FastaSeqFetcher) can differ from the RefSeq transcript,
  eg the GRCh37 reference has the ACTN3 R577X stop allele, which makes every p. of the
  transcript "p.?"
* report exceptions and conversion problems cdot knows about but can't apply

Problems are returned as HGVSFix, as with fix_hgvs().
"""

import copy
from typing import Optional

from bioutils.sequences import TranslationTable, aa3_to_aa1, translate_cds
from hgvs.edit import AAExt, AAFs, AASub
from hgvs.exceptions import HGVSUnsupportedOperationError

from cdot.hgvs.clean import HGVSFix, HGVSFixCode, HGVSFixSeverity

# NCBI genetic code number -> (biocommons translation table, name)
_TRANSLATION_TABLES = {
    1: (TranslationTable.standard, "standard"),
    2: (TranslationTable.vertebrate_mitochondrial, "vertebrate mitochondrial"),
}

# RefSeq CDS exceptions fix_c_to_p deals with (or checks by codon) so doesn't report as not applied
_EXCEPTION_RIBOSOMAL_SLIPPAGE = "ribosomal slippage"
_START_CODON_EXCEPTIONS = {"alternative start codon", "translation initiation by tRNA-Leu at CUG codon"}


def fix_c_to_p(variant_mapper, var_c, raise_on_errors: bool = False, **kwargs) -> tuple:
    """
    Convert a c. SequenceVariant to p. with a biocommons VariantMapper or AssemblyMapper,
    using cdot's translation data (from the mapper's data provider) to choose the genetic code,
    fix what it can and report what it can't.

    kwargs are passed through to variant_mapper.c_to_p() (eg pro_ac/alt_ac for a VariantMapper).
    If translation_table is passed, it is used as is rather than chosen from the data.

    Returns:
        (var_p, fixes)
        var_p is None if there was an ERROR-level fix, as the p. would be wrong.
        With raise_on_errors=True, raises HGVSUnsupportedOperationError on the first ERROR instead.

    Transcripts without translation data (eg data older than schema 0.2.35, or a data
    provider that isn't cdot's) are converted as biocommons does, with no fixes.

    Example::

        var_p, fixes = fix_c_to_p(am, hp.parse("NM_002085.5:c.217T>C"))
        # var_p = NP_002076.2:p.(Sec73Arg)
        # fixes = [HGVSFix(WARNING, USED_TRANSLATION_TABLE, ...)]
    """
    transcript = None
    if get_transcript := getattr(variant_mapper.hdp, "_get_transcript", None):
        transcript = get_transcript(var_c.ac)

    fixes = []
    sec_codons = set()
    translation_table = kwargs.get("translation_table")
    if transcript is not None:
        translation = transcript.get("translation") or {}
        sec_codons = set((translation.get("transl_except") or {}).get("Sec", []))
        if translation_table is None:
            translation_table, table_fixes = _get_translation_table(var_c.ac, translation, sec_codons)
            fixes.extend(table_fixes)
        fixes.extend(_slippage_fixes(var_c.ac, transcript, translation))
        fixes.extend(_build_warning_fixes(var_c.ac, transcript))
        fixes.extend(_exception_fixes(var_c.ac, translation))
        if translation_table is not None and not _errors(fixes):
            fixes.extend(_codon_fixes(variant_mapper.hdp, var_c.ac, transcript, translation, translation_table))

    if errors := _errors(fixes):
        if raise_on_errors:
            raise HGVSUnsupportedOperationError(errors[0].message)
        return None, fixes

    if translation_table is not None:
        kwargs["translation_table"] = translation_table
    var_p = variant_mapper.c_to_p(var_c, **kwargs)

    if translation_table == TranslationTable.selenocysteine:
        var_p, sec_fixes = _fix_selenoprotein_p(var_p, sec_codons)
        fixes.extend(sec_fixes)
    return var_p, fixes


def _errors(fixes: list[HGVSFix]) -> list[HGVSFix]:
    return [f for f in fixes if f.severity == HGVSFixSeverity.ERROR]


def _get_translation_table(tx_ac: str, translation, sec_codons: set) -> tuple:
    fixes = []
    transl_table = translation.get("transl_table") or 1
    translation_table, table_name = _TRANSLATION_TABLES.get(transl_table, (None, None))
    if translation_table is None:
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.ERROR,
            code=HGVSFixCode.UNSUPPORTED_TRANSLATION_TABLE,
            message=f"{tx_ac} uses genetic code {transl_table}, which biocommons HGVS can't translate",
        ))
    elif transl_table != 1:
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.USED_TRANSLATION_TABLE,
            message=f"Translated {tx_ac} with the {table_name} genetic code ({transl_table})",
        ))
    elif sec_codons:
        translation_table = TranslationTable.selenocysteine
        codons = ", ".join(str(c) for c in sorted(sec_codons))
        fixes.append(HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.USED_TRANSLATION_TABLE,
            message=f"{tx_ac} is a selenoprotein (Sec at codon {codons}), translated with UGA as Sec",
        ))
    return translation_table, fixes


def _slippage_fixes(tx_ac: str, transcript, translation) -> list[HGVSFix]:
    has_slippage = bool(translation.get("ribosomal_slippage")) \
        or _EXCEPTION_RIBOSOMAL_SLIPPAGE in (translation.get("exceptions") or []) \
        or any((build_data.get("warnings") or {}).get("ribosomal_slippage_unplaced")
               for build_data in transcript.get("genome_builds", {}).values())
    if not has_slippage:
        return []
    return [HGVSFix(
        severity=HGVSFixSeverity.ERROR,
        code=HGVSFixCode.RIBOSOMAL_SLIPPAGE_UNSUPPORTED,
        message=f"{tx_ac} has a programmed ribosomal frameshift, which biocommons HGVS can't translate",
    )]


def _build_warning_fixes(tx_ac: str, transcript) -> list[HGVSFix]:
    """ Problems cdot hit converting the transcript (ribosomal_slippage_unplaced is an ERROR elsewhere) """
    fixes = []
    for genome_build, build_data in transcript.get("genome_builds", {}).items():
        warnings = build_data.get("warnings") or {}
        if amino_acids := warnings.get("transl_except_unplaced"):
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=HGVSFixCode.TRANSLATION_DATA_INCOMPLETE,
                message=f"{tx_ac} has translation exceptions for {', '.join(amino_acids)} on {genome_build} "
                        "that cdot couldn't place, so couldn't apply",
            ))
        if codons := warnings.get("codons_unplaced"):
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=HGVSFixCode.TRANSLATION_DATA_INCOMPLETE,
                message=f"{tx_ac} is missing {' and '.join(codons)} on {genome_build}",
            ))
    return fixes


def _exception_fixes(tx_ac: str, translation) -> list[HGVSFix]:
    """ RefSeq CDS exceptions that cdot doesn't interpret """
    handled = {_EXCEPTION_RIBOSOMAL_SLIPPAGE}
    if 1 in _transl_except_codons(translation.get("transl_except")):
        handled |= _START_CODON_EXCEPTIONS  # checked by codon
    fixes = []
    for exception in translation.get("exceptions") or []:
        if exception not in handled:
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=HGVSFixCode.TRANSLATION_EXCEPTION_NOT_APPLIED,
                message=f"{tx_ac} has RefSeq CDS exception '{exception}', the p. may not match the RefSeq protein",
            ))
    return fixes


def _transl_except_codons(transl_except) -> dict:
    """ {codon: amino_acid} """
    return {codon: amino_acid for amino_acid, codons in (transl_except or {}).items() for codon in codons}


def _codon_fixes(hdp, tx_ac: str, transcript, translation, translation_table) -> list[HGVSFix]:
    """ Check the codons RefSeq gives a translation exception for against the transcript sequence
        biocommons will translate. Genome mismatches of every build are checked, as the sequence
        can be from any build's genome (or the RefSeq transcript, when they all match) """
    codon_amino_acids = _transl_except_codons(translation.get("transl_except"))
    genome_mismatch_codons = set()
    for build_data in transcript.get("genome_builds", {}).values():
        genome_mismatch = build_data.get("genome_mismatch") or {}
        mismatch_codon_amino_acids = _transl_except_codons(genome_mismatch.get("transl_except"))
        codon_amino_acids.update(mismatch_codon_amino_acids)
        genome_mismatch_codons.update(mismatch_codon_amino_acids)
    codon_amino_acids = {c: aa for c, aa in codon_amino_acids.items() if aa != "TERM"}  # biocommons pads mito
    if not codon_amino_acids:
        return []

    tx_info = hdp.get_tx_identity_info(tx_ac)
    cds = hdp.get_seq(tx_ac)[tx_info["cds_start_i"]:tx_info["cds_end_i"]]
    last_codon = len(cds) // 3
    fixes = []
    for codon, amino_acid in sorted(codon_amino_acids.items()):
        codon_seq = cds[(codon - 1) * 3:codon * 3]
        if len(codon_seq) != 3:
            continue
        translated_aa = translate_cds(codon_seq, translation_table=translation_table)
        if amino_acid == "Other":
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=HGVSFixCode.TRANSLATION_EXCEPTION_NOT_APPLIED,
                message=f"{tx_ac} codon {codon} ({codon_seq}) is a non-standard amino acid in the RefSeq protein "
                        f"(eg stop codon readthrough), but is translated as {_aa3(translated_aa)}",
            ))
            continue
        expected_aa = aa3_to_aa1(amino_acid)
        if translated_aa == expected_aa:
            continue
        if translated_aa == "*" and codon != last_codon:
            message = f"{tx_ac} codon {codon} is {codon_seq}, translated as a stop, but is {amino_acid} in the " \
                      "RefSeq protein, so every p. would be wrong"
            if codon in genome_mismatch_codons:
                message += ". The transcript sequence is probably from a genome that differs from the " \
                           "transcript here, use the RefSeq transcript sequence (eg SeqRepo)"
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.ERROR,
                code=HGVSFixCode.REFERENCE_CODON_MISMATCH,
                message=message,
            ))
        else:
            fixes.append(HGVSFix(
                severity=HGVSFixSeverity.WARNING,
                code=HGVSFixCode.REFERENCE_CODON_MISMATCH,
                message=f"{tx_ac} codon {codon} is {codon_seq} ({_aa3(translated_aa)}) in the transcript sequence "
                        f"used, but {amino_acid} in the RefSeq protein",
            ))
    return fixes


def _aa3(aa1: str) -> str:
    return "Ter" if aa1 == "*" else _AA1_TO_AA3.get(aa1, aa1)


_AA1_TO_AA3 = {aa3_to_aa1(aa3): aa3 for aa3 in (
    "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met",
    "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "Sec",
)}


def _fix_selenoprotein_p(var_p, sec_codons: set) -> tuple:
    """ Translating with the selenocysteine code reads every UGA as Sec, but only the known codons
        are (they need a SECIS element), so a new UGA is a stop """
    posedit = var_p.posedit
    if posedit is None or posedit.edit is None:
        return var_p, []
    edit = posedit.edit
    start = posedit.pos.start.base if posedit.pos else None
    if isinstance(edit, AASub) and edit.alt == "U" and start not in sec_codons:
        fixed_p = copy.deepcopy(var_p)
        fixed_p.posedit.edit = AASub(ref=edit.ref, alt="*", uncertain=edit.uncertain, init_met=edit.init_met)
        return fixed_p, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.NEW_UGA_READ_AS_STOP,
            message=f"The variant makes a new UGA at codon {start} of a selenoprotein, which is a stop not Sec",
            original=str(var_p),
            fixed=str(fixed_p),
        )]
    if isinstance(edit, (AAFs, AAExt)) or "U" in (getattr(edit, "alt", None) or ""):
        return var_p, [HGVSFix(
            severity=HGVSFixSeverity.WARNING,
            code=HGVSFixCode.UGA_MAY_BE_READ_AS_SEC,
            message="The variant changes the sequence translated in a selenoprotein, where UGA is read as Sec, "
                    "so a new UGA stop may have been read as Sec",
        )]
    return var_p, []
