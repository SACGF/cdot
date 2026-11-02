"""
Tests for cdot.hgvs.translation.fix_c_to_p - c. to p. using cdot's translation data (#131, #76).

Transcript data is from the translation test GFFs (see test_gff_parsers). The transcript
sequences are spliced from the genome, as FastaSeqFetcher would, so ACTN3 on GRCh37 has the
reference R577X stop allele.
"""
import json
import os

import hgvs.parser
import pytest
from bioutils.sequences import TranslationTable
from hgvs.assemblymapper import AssemblyMapper
from hgvs.exceptions import HGVSUnsupportedOperationError

from cdot.hgvs.clean import HGVSFixCode, HGVSFixSeverity
from cdot.hgvs.dataproviders.json_data_provider import AbstractJSONDataProvider, JSONDataProvider
from cdot.hgvs.translation import (
    _build_warning_fixes,
    _codon_fixes,
    _exception_fixes,
    _fix_selenoprotein_p,
    _get_translation_table,
    _slippage_fixes,
    fix_c_to_p,
)
from tests.mock_seqfetcher import MockSeqFetcher

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "test_data")
SEQUENCES_JSON = os.path.join(TEST_DATA_DIR, "translation_transcript_sequences.json")

C = HGVSFixCode
W = HGVSFixSeverity.WARNING
E = HGVSFixSeverity.ERROR

hp = hgvs.parser.Parser()


def _assembly_mapper(genome_build, seqfetcher=None):
    json_file = os.path.join(TEST_DATA_DIR, f"cdot.refseq.translation.{genome_build.lower()}.json")
    hdp = JSONDataProvider([json_file], seqfetcher=seqfetcher or MockSeqFetcher(SEQUENCES_JSON))
    return AssemblyMapper(hdp, assembly_name=genome_build, alt_aln_method="splign", replace_reference=True)


@pytest.fixture(scope="module")
def am_grch38():
    return _assembly_mapper("GRCh38")


@pytest.fixture(scope="module")
def am_grch37():
    return _assembly_mapper("GRCh37")


def _codes(fixes):
    return [(f.severity, f.code) for f in fixes]


class TestSelenoprotein:
    def test_biocommons_alone_reads_sec_as_stop(self, am_grch38):
        """ Why the wrapper is needed: without the protein sequence, biocommons translates GPX4's Sec as a stop """
        assert str(am_grch38.c_to_p(hp.parse("NM_002085.5:c.100G>A"))) == "NP_002076.2:p.?"

    def test_missense(self, am_grch38):
        var_p, fixes = fix_c_to_p(am_grch38, hp.parse("NM_002085.5:c.100G>A"))
        assert str(var_p) == "NP_002076.2:p.(Asp34Asn)"
        assert _codes(fixes) == [(W, C.USED_TRANSLATION_TABLE)]

    def test_sec_codon(self, am_grch38):
        var_p, _ = fix_c_to_p(am_grch38, hp.parse("NM_002085.5:c.217T>C"))
        assert str(var_p) == "NP_002076.2:p.(Sec73Arg)"

    def test_new_uga_is_stop(self, am_grch38):
        var_p, fixes = fix_c_to_p(am_grch38, hp.parse("NM_002085.5:c.105G>A"))
        assert str(var_p) == "NP_002076.2:p.(Trp35Ter)"
        assert _codes(fixes) == [(W, C.USED_TRANSLATION_TABLE), (W, C.NEW_UGA_READ_AS_STOP)]
        assert fixes[-1].original == "NP_002076.2:p.(Trp35Sec)"
        assert fixes[-1].fixed == "NP_002076.2:p.(Trp35Ter)"

    def test_frameshift_warns(self, am_grch38):
        var_p, fixes = fix_c_to_p(am_grch38, hp.parse("NM_002085.5:c.105del"))
        assert "fs" in str(var_p)
        assert (W, C.UGA_MAY_BE_READ_AS_SEC) in _codes(fixes)

    def test_explicit_translation_table_is_used(self, am_grch38):
        """ Forcing the standard code reads Sec as a stop, which the codon check catches """
        var_c = hp.parse("NM_002085.5:c.100G>A")
        var_p, fixes = fix_c_to_p(am_grch38, var_c, translation_table=TranslationTable.standard)
        assert var_p is var_c  # the input, unchanged, as fix_hgvs does on an error
        assert _codes(fixes) == [(E, C.REFERENCE_CODON_MISMATCH)]
        assert "SeqRepo" not in fixes[0].message


def test_mitochondrial(am_grch38):
    var_p, fixes = fix_c_to_p(am_grch38, hp.parse("fake-rna-COX1:c.100A>G"))
    assert str(var_p) == "YP_003024028.1:p.(Ser34Gly)"
    assert _codes(fixes) == [(W, C.USED_TRANSLATION_TABLE)]


class TestRibosomalSlippage:
    def test_error(self, am_grch38):
        var_c = hp.parse("NM_015068.3:c.100G>A")
        var_p, fixes = fix_c_to_p(am_grch38, var_c)
        assert var_p is var_c
        assert _codes(fixes) == [(E, C.RIBOSOMAL_SLIPPAGE_UNSUPPORTED)]

    def test_raise_on_errors(self, am_grch38):
        with pytest.raises(HGVSUnsupportedOperationError):
            fix_c_to_p(am_grch38, hp.parse("NM_015068.3:c.100G>A"), raise_on_errors=True)


class TestGenomeMismatch:
    def test_genome_stop_codon_is_error(self, am_grch37):
        """ GRCh37 has the ACTN3 R577X stop allele, so the genome derived sequence has a premature stop """
        var_c = hp.parse("NM_001104.4:c.100G>A")
        var_p, fixes = fix_c_to_p(am_grch37, var_c)
        assert var_p is var_c
        assert _codes(fixes) == [(E, C.REFERENCE_CODON_MISMATCH)]
        assert "codon 577" in fixes[0].message
        assert "SeqRepo" in fixes[0].message

    def test_refseq_sequence_is_fine(self, tmp_path):
        """ With the RefSeq transcript sequence (c.1729 is C, R577), there's nothing to report """
        with open(SEQUENCES_JSON) as f:
            sequences = json.load(f)
        am = _assembly_mapper("GRCh37")
        cds_start_i = am.hdp.get_tx_identity_info("NM_001104.4")["cds_start_i"]
        seq = sequences["NM_001104.4"]
        i = cds_start_i + 1728
        assert seq[i] == "T"
        sequences["NM_001104.4"] = seq[:i] + "C" + seq[i + 1:]
        refseq_sequences = tmp_path / "sequences.json"
        refseq_sequences.write_text(json.dumps(sequences))

        am = _assembly_mapper("GRCh37", seqfetcher=MockSeqFetcher(str(refseq_sequences)))
        var_p, fixes = fix_c_to_p(am, hp.parse("NM_001104.4:c.100G>A"))
        assert str(var_p) == "NP_001095.2:p.(Asp34Asn)"
        assert fixes == []


class TestTranscriptRecordAccess:
    def test_provider_get_transcript(self, am_grch38):
        """ The public accessor fix_c_to_p reads the record through """
        record = am_grch38.hdp.get_transcript("NM_002085.5")
        assert record["id"] == "NM_002085.5"
        assert record["translation"]["transl_except"] == {"Sec": [73]}
        assert am_grch38.hdp.get_transcript("NM_999999.9") is None

    def test_fix_c_to_p_uses_public_accessor(self, am_grch38, monkeypatch):
        calls = []
        real_get_transcript = AbstractJSONDataProvider.get_transcript

        def recording_get_transcript(self, tx_ac):
            calls.append(tx_ac)
            return real_get_transcript(self, tx_ac)

        monkeypatch.setattr(AbstractJSONDataProvider, "get_transcript", recording_get_transcript)
        fix_c_to_p(am_grch38, hp.parse("NM_002085.5:c.100G>A"))
        assert calls == ["NM_002085.5"]

    def test_provider_without_accessor_converts_as_biocommons(self, am_grch38, monkeypatch):
        """ A data provider that isn't cdot's has no get_transcript, so there is no translation data """
        monkeypatch.delattr(AbstractJSONDataProvider, "get_transcript")
        var_p, fixes = fix_c_to_p(am_grch38, hp.parse("NM_002085.5:c.100G>A"))
        assert str(var_p) == "NP_002076.2:p.?"
        assert fixes == []


def test_exception_fixes():
    translation = {"exceptions": ["unclassified translation discrepancy"]}
    assert _codes(_exception_fixes("NM_1.1", translation)) == [(W, C.TRANSLATION_EXCEPTION_NOT_APPLIED)]
    # A start codon exception is checked by codon when there's a codon 1 transl_except
    translation = {"exceptions": ["alternative start codon"], "transl_except": {"Met": [1]}}
    assert _exception_fixes("NM_1.1", translation) == []
    translation = {"exceptions": ["alternative start codon"]}
    assert _codes(_exception_fixes("NM_1.1", translation)) == [(W, C.TRANSLATION_EXCEPTION_NOT_APPLIED)]


def test_build_warning_fixes():
    transcript = {"genome_builds": {
        "GRCh38": {"warnings": {"transl_except_unplaced": ["Sec"], "codons_unplaced": ["start_codon"]}},
    }}
    fixes = _build_warning_fixes("ENST1.1", transcript)
    assert _codes(fixes) == [(W, C.TRANSLATION_DATA_INCOMPLETE)] * 2
    assert "GRCh38" in fixes[0].message


class TestGetTranslationTable:
    def test_default_is_standard(self):
        assert _get_translation_table("NM_000059.4", {}, set()) == (TranslationTable.standard, [])
        assert _get_translation_table("NM_000059.4", {"transl_table": 1}, set()) == (TranslationTable.standard, [])

    def test_mitochondrial(self):
        table, fixes = _get_translation_table("fake-rna-COX1", {"transl_table": 2}, set())
        assert table == TranslationTable.vertebrate_mitochondrial
        assert _codes(fixes) == [(W, C.USED_TRANSLATION_TABLE)]
        assert "vertebrate mitochondrial genetic code (2)" in fixes[0].message

    def test_unsupported_table_is_error(self):
        table, fixes = _get_translation_table("NM_000059.4", {"transl_table": 3}, set())
        assert table is None
        assert _codes(fixes) == [(E, C.UNSUPPORTED_TRANSLATION_TABLE)]
        assert "NM_000059.4 uses genetic code 3" in fixes[0].message

    def test_selenoprotein(self):
        table, fixes = _get_translation_table("NM_002085.5", {}, {73})
        assert table == TranslationTable.selenocysteine
        assert _codes(fixes) == [(W, C.USED_TRANSLATION_TABLE)]
        assert "Sec at codon 73" in fixes[0].message

    def test_selenoprotein_codons_sorted_in_message(self):
        _table, fixes = _get_translation_table("NM_002085.5", {}, {120, 73})
        assert "Sec at codon 73, 120" in fixes[0].message

    def test_non_standard_table_takes_precedence_over_sec_codons(self):
        table, fixes = _get_translation_table("NM_000059.4", {"transl_table": 2}, {73})
        assert table == TranslationTable.vertebrate_mitochondrial
        assert _codes(fixes) == [(W, C.USED_TRANSLATION_TABLE)]


class TestUnsupportedTranslationTable:
    """ A genetic code biocommons has no table for: no p. is produced """

    @pytest.fixture
    def am_table_3(self, am_grch38, monkeypatch):
        record = {"id": "NM_002085.5", "translation": {"transl_table": 3}, "genome_builds": {}}
        monkeypatch.setattr(am_grch38.hdp, "get_transcript", lambda tx_ac: record)
        return am_grch38

    def test_returns_input_with_error(self, am_table_3):
        var_c = hp.parse("NM_002085.5:c.100G>A")
        var_p, fixes = fix_c_to_p(am_table_3, var_c)
        assert var_p is var_c
        assert _codes(fixes) == [(E, C.UNSUPPORTED_TRANSLATION_TABLE)]

    def test_raise_on_errors(self, am_table_3):
        with pytest.raises(HGVSUnsupportedOperationError, match="genetic code 3"):
            fix_c_to_p(am_table_3, hp.parse("NM_002085.5:c.100G>A"), raise_on_errors=True)

    def test_explicit_table_overrides(self, am_table_3):
        """ A caller supplied translation_table is used as is, so the data's code is not looked at """
        var_p, fixes = fix_c_to_p(am_table_3, hp.parse("NM_002085.5:c.100G>A"),
                                  translation_table=TranslationTable.standard)
        assert str(var_p) == "NP_002076.2:p.?"  # standard code reads the Sec codon as a stop
        assert fixes == []


class TestSlippageFixes:
    def test_no_slippage(self):
        assert _slippage_fixes("NM_000059.4", {"genome_builds": {"GRCh38": {}}}, {}) == []
        assert _slippage_fixes("NM_000059.4", {"genome_builds": {"GRCh38": {"warnings": {}}}},
                               {"exceptions": ["unclassified translation discrepancy"]}) == []

    def test_slippage_data(self):
        translation = {"ribosomal_slippage": [{"cds_position": 957, "shift": -1}]}
        fixes = _slippage_fixes("NM_015068.3", {"genome_builds": {}}, translation)
        assert _codes(fixes) == [(E, C.RIBOSOMAL_SLIPPAGE_UNSUPPORTED)]
        assert "NM_015068.3 has a programmed ribosomal frameshift" in fixes[0].message

    def test_slippage_exception_only(self):
        fixes = _slippage_fixes("NM_015068.3", {"genome_builds": {}}, {"exceptions": ["ribosomal slippage"]})
        assert _codes(fixes) == [(E, C.RIBOSOMAL_SLIPPAGE_UNSUPPORTED)]

    def test_slippage_unplaced_build_warning(self):
        """ cdot knew about the slippage but couldn't place it: still an error, not a build warning """
        transcript = {"genome_builds": {"GRCh38": {"warnings": {"ribosomal_slippage_unplaced": True}}}}
        fixes = _slippage_fixes("NM_015068.3", transcript, {})
        assert _codes(fixes) == [(E, C.RIBOSOMAL_SLIPPAGE_UNSUPPORTED)]
        assert _build_warning_fixes("NM_015068.3", transcript) == []

    def test_slippage_unplaced_in_one_build_of_several(self):
        transcript = {"genome_builds": {
            "GRCh37": {"warnings": None},
            "GRCh38": {"warnings": {"ribosomal_slippage_unplaced": True}},
        }}
        assert _codes(_slippage_fixes("NM_015068.3", transcript, {})) == [(E, C.RIBOSOMAL_SLIPPAGE_UNSUPPORTED)]

    def test_slippage_reported_once(self):
        transcript = {"genome_builds": {"GRCh38": {"warnings": {"ribosomal_slippage_unplaced": True}}}}
        translation = {"ribosomal_slippage": [{"cds_position": 957, "shift": -1}],
                       "exceptions": ["ribosomal slippage"]}
        assert len(_slippage_fixes("NM_015068.3", transcript, translation)) == 1


class _CdsProvider:
    """ Just enough of a data provider for _codon_fixes: the sequence is the CDS """

    def __init__(self, cds):
        self.cds = cds
        self.calls = 0

    def get_tx_identity_info(self, tx_ac):
        self.calls += 1
        return {"cds_start_i": 0, "cds_end_i": len(self.cds)}

    def get_seq(self, tx_ac):
        return self.cds


class TestCodonFixes:
    # Met, TGA (a stop in the standard code), Gly, Trp, stop: 5 codons
    CDS = "ATG" "TGA" "GGG" "TGG" "TAA"

    def _fixes(self, transl_except, genome_mismatch=None, translation_table=TranslationTable.standard):
        hdp = _CdsProvider(self.CDS)
        build = {"genome_mismatch": {"transl_except": genome_mismatch}} if genome_mismatch else {}
        transcript = {"genome_builds": {"GRCh38": build}}
        return hdp, _codon_fixes(hdp, "NM_000059.4", transcript, {"transl_except": transl_except}, translation_table)

    def test_no_exceptions_no_sequence_fetch(self):
        hdp, fixes = self._fixes(None)
        assert fixes == []
        assert hdp.calls == 0

    def test_term_exceptions_are_ignored(self):
        """ A TERM transl_except is biocommons padding a mito stop, there's no amino acid to compare """
        hdp, fixes = self._fixes({"TERM": [5]})
        assert fixes == []
        assert hdp.calls == 0
        hdp, fixes = self._fixes({"TERM": [5], "Gly": [3]})
        assert fixes == []
        assert hdp.calls == 1

    def test_matching_codons_no_fix(self):
        _hdp, fixes = self._fixes({"Met": [1], "Gly": [3], "Trp": [4]})
        assert fixes == []

    def test_other_amino_acid_is_warning(self):
        """ RefSeq's "Other" (eg a readthrough) can't be checked, but the translation is reported """
        _hdp, fixes = self._fixes({"Other": [2]})
        assert _codes(fixes) == [(W, C.TRANSLATION_EXCEPTION_NOT_APPLIED)]
        assert "codon 2 (TGA) is a non-standard amino acid" in fixes[0].message
        assert "translated as Ter" in fixes[0].message

    def test_non_stop_mismatch_is_warning(self):
        _hdp, fixes = self._fixes({"Met": [3]})
        assert _codes(fixes) == [(W, C.REFERENCE_CODON_MISMATCH)]
        assert "codon 3 is GGG (Gly)" in fixes[0].message
        assert "but Met in the RefSeq protein" in fixes[0].message

    def test_internal_stop_mismatch_is_error(self):
        _hdp, fixes = self._fixes({"Trp": [2]})
        assert _codes(fixes) == [(E, C.REFERENCE_CODON_MISMATCH)]
        assert "codon 2 is TGA, translated as a stop, but is Trp" in fixes[0].message
        assert "every p. would be wrong" in fixes[0].message
        assert "SeqRepo" not in fixes[0].message

    def test_stop_at_last_codon_is_warning(self):
        """ A stop where the CDS ends doesn't truncate the protein """
        _hdp, fixes = self._fixes({"Trp": [5]})
        assert _codes(fixes) == [(W, C.REFERENCE_CODON_MISMATCH)]
        assert "codon 5 is TAA (Ter)" in fixes[0].message

    def test_genome_mismatch_codon_blames_the_sequence(self):
        _hdp, fixes = self._fixes(None, genome_mismatch={"Trp": [2]})
        assert _codes(fixes) == [(E, C.REFERENCE_CODON_MISMATCH)]
        assert "SeqRepo" in fixes[0].message

    def test_genome_mismatch_merged_with_transcript_exceptions(self):
        _hdp, fixes = self._fixes({"Met": [3]}, genome_mismatch={"Trp": [2]})
        assert _codes(fixes) == [(E, C.REFERENCE_CODON_MISMATCH), (W, C.REFERENCE_CODON_MISMATCH)]
        assert "codon 2" in fixes[0].message and "codon 3" in fixes[1].message

    def test_codon_beyond_cds_skipped(self):
        _hdp, fixes = self._fixes({"Gly": [9]})
        assert fixes == []

    def test_mitochondrial_table_reads_tga_as_trp(self):
        _hdp, fixes = self._fixes({"Trp": [2]}, translation_table=TranslationTable.vertebrate_mitochondrial)
        assert fixes == []


class TestFixSelenoproteinP:
    def _fix(self, hgvs_p, sec_codons=frozenset()):
        return _fix_selenoprotein_p(hp.parse(hgvs_p), set(sec_codons))

    def test_new_sec_is_stop(self):
        var_p, fixes = self._fix("NP_002076.2:p.(Trp35Sec)")
        assert str(var_p) == "NP_002076.2:p.(Trp35Ter)"
        assert _codes(fixes) == [(W, C.NEW_UGA_READ_AS_STOP)]
        assert fixes[0].original == "NP_002076.2:p.(Trp35Sec)"
        assert fixes[0].fixed == "NP_002076.2:p.(Trp35Ter)"

    def test_sec_at_known_codon_not_rewritten(self):
        """ At a codon RefSeq says is Sec the UGA really is Sec, so it isn't turned into a stop """
        var_p, fixes = self._fix("NP_002076.2:p.(Arg73Sec)", sec_codons={73})
        assert str(var_p) == "NP_002076.2:p.(Arg73Sec)"
        assert _codes(fixes) == [(W, C.UGA_MAY_BE_READ_AS_SEC)]

    def test_extension_may_read_uga_as_sec(self):
        """ The extension past the stop could contain a UGA biocommons translated as Sec """
        var_p, fixes = self._fix("NP_002076.2:p.(Ter200GlnextTer10)")
        assert str(var_p) == "NP_002076.2:p.(Ter200GlnextTer10)"
        assert _codes(fixes) == [(W, C.UGA_MAY_BE_READ_AS_SEC)]

    def test_frameshift_may_read_uga_as_sec(self):
        _var_p, fixes = self._fix("NP_002076.2:p.(Asp34fs)")
        assert _codes(fixes) == [(W, C.UGA_MAY_BE_READ_AS_SEC)]

    def test_delins_with_sec_may_read_uga_as_sec(self):
        _var_p, fixes = self._fix("NP_002076.2:p.(Asp34_Trp35delinsSecGly)")
        assert _codes(fixes) == [(W, C.UGA_MAY_BE_READ_AS_SEC)]

    def test_delins_without_sec_no_fix(self):
        _var_p, fixes = self._fix("NP_002076.2:p.(Asp34_Trp35delinsGly)")
        assert fixes == []

    @pytest.mark.parametrize("hgvs_p", [
        "NP_002076.2:p.(Asp34Asn)",
        "NP_002076.2:p.(Met1?)",
        "NP_002076.2:p.?",
        "NP_002076.2:p.=",
    ])
    def test_nothing_to_do(self, hgvs_p):
        var_p, fixes = self._fix(hgvs_p)
        assert str(var_p) == hgvs_p
        assert fixes == []
