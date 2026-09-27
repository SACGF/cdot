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
from cdot.hgvs.dataproviders.json_data_provider import JSONDataProvider
from cdot.hgvs.translation import _build_warning_fixes, _exception_fixes, fix_c_to_p
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
        var_p, fixes = fix_c_to_p(am_grch38, hp.parse("NM_002085.5:c.100G>A"),
                              translation_table=TranslationTable.standard)
        assert var_p is None
        assert _codes(fixes) == [(E, C.REFERENCE_CODON_MISMATCH)]
        assert "SeqRepo" not in fixes[0].message


def test_mitochondrial(am_grch38):
    var_p, fixes = fix_c_to_p(am_grch38, hp.parse("fake-rna-COX1:c.100A>G"))
    assert str(var_p) == "YP_003024028.1:p.(Ser34Gly)"
    assert _codes(fixes) == [(W, C.USED_TRANSLATION_TABLE)]


class TestRibosomalSlippage:
    def test_error(self, am_grch38):
        var_p, fixes = fix_c_to_p(am_grch38, hp.parse("NM_015068.3:c.100G>A"))
        assert var_p is None
        assert _codes(fixes) == [(E, C.RIBOSOMAL_SLIPPAGE_UNSUPPORTED)]

    def test_raise_on_errors(self, am_grch38):
        with pytest.raises(HGVSUnsupportedOperationError):
            fix_c_to_p(am_grch38, hp.parse("NM_015068.3:c.100G>A"), raise_on_errors=True)


class TestGenomeMismatch:
    def test_genome_stop_codon_is_error(self, am_grch37):
        """ GRCh37 has the ACTN3 R577X stop allele, so the genome derived sequence has a premature stop """
        var_p, fixes = fix_c_to_p(am_grch37, hp.parse("NM_001104.4:c.100G>A"))
        assert var_p is None
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
