"""
Tests for FastaSeqFetcher warning/refusing to build a transcript that differs from the genome (#115)

Uses a tiny synthetic genome and transcripts. The genome_mismatch values are from the public examples
in #115: EGR2 NM_000399.3 (1 indel) and PTEN NM_000314.4 (1 substitution, 1 indel)
"""
import json
import warnings

import pysam
import pytest
from hgvs.exceptions import HGVSDataNotAvailableError

from cdot.hgvs.dataproviders import ChainedSeqFetcher, JSONDataProvider
from cdot.hgvs.dataproviders.fasta_seqfetcher import FastaSeqFetcher, ExonsFromGenomeFastaSeqFetcher, \
    GenomeMismatchError, GenomeMismatchPolicy, GenomeMismatchWarning

CONTIG = "NC_000011.9"  # GRCh37 chr11
GENOME = "ACGTTGCAAC" * 10
EXONS = [[10, 20, 0, 1, 10, None], [30, 40, 1, 11, 20, None]]
TRANSCRIPT_SEQ = GENOME[10:20] + GENOME[30:40]

GENOME_MISMATCH = {
    "NM_000399.3": {"exceptions": ["annotated by transcript or proteomic data"],
                    "transcript": {"non_frameshifting_indels": 1},
                    "alignment": {"num_mismatch": 0, "gap_count": 1, "pct_identity_gap": 99.9664,
                                  "pct_coverage": 100}},
    "NM_000314.4": {"exceptions": ["annotated by transcript or proteomic data"],
                    "transcript": {"substitutions": 1, "non_frameshifting_indels": 1}},
}
MATCHES_GENOME = "NM_000059.4"


def _transcript(accession, genome_mismatch=None):
    build = {"contig": CONTIG, "strand": "+", "exons": EXONS}
    if genome_mismatch:
        build["genome_mismatch"] = genome_mismatch
    return {"id": accession, "gene_name": "TEST", "biotype": ["ncRNA"], "genome_builds": {"GRCh37": build}}


@pytest.fixture(scope="module")
def test_files(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("fasta_seqfetcher")
    fasta_filename = str(tmp_path / "genome.fa")
    with open(fasta_filename, "w") as f:
        f.write(f">{CONTIG}\n{GENOME}\n")
    pysam.faidx(fasta_filename)

    transcripts = {ac: _transcript(ac, gm) for ac, gm in GENOME_MISMATCH.items()}
    transcripts[MATCHES_GENOME] = _transcript(MATCHES_GENOME)
    json_filename = str(tmp_path / "cdot.json")
    with open(json_filename, "w") as f:
        json.dump({"cdot_version": "0.2.36", "genome_builds": ["GRCh37"],
                   "transcripts": transcripts, "genes": {}}, f)
    return fasta_filename, json_filename


def _hdp(json_filename, seqfetcher):
    return JSONDataProvider([json_filename], seqfetcher=seqfetcher)


class StaticSeqFetcher:
    source = "static"

    def fetch_seq(self, ac, start_i=None, end_i=None):
        return "REAL"


@pytest.mark.parametrize("seqfetcher_class", [FastaSeqFetcher, ExonsFromGenomeFastaSeqFetcher])
def test_warns_by_default(test_files, seqfetcher_class):
    fasta_filename, json_filename = test_files
    hdp = _hdp(json_filename, seqfetcher_class(fasta_filename))
    with pytest.warns(GenomeMismatchWarning, match=r"NM_000314.4 differs from the genome \(NC_000011.9\)"):
        assert hdp.seqfetcher.fetch_seq("NM_000314.4") == TRANSCRIPT_SEQ


def test_warning_describes_mismatch(test_files):
    fasta_filename, json_filename = test_files
    hdp = _hdp(json_filename, FastaSeqFetcher(fasta_filename))
    with pytest.warns(GenomeMismatchWarning) as record:
        hdp.seqfetcher.fetch_seq("NM_000399.3")
    message = str(record[0].message)
    assert "non_frameshifting_indels" in message
    assert "annotated by transcript or proteomic data" in message


def test_warns_once_per_transcript(test_files):
    fasta_filename, json_filename = test_files
    hdp = _hdp(json_filename, FastaSeqFetcher(fasta_filename, cache=False))
    with pytest.warns(GenomeMismatchWarning) as record:
        hdp.seqfetcher.fetch_seq("NM_000399.3")
        hdp.seqfetcher.fetch_seq("NM_000399.3", 0, 5)
    assert len(record) == 1


def test_no_warning_when_transcript_matches_genome(test_files):
    fasta_filename, json_filename = test_files
    hdp = _hdp(json_filename, FastaSeqFetcher(fasta_filename))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert hdp.seqfetcher.fetch_seq(MATCHES_GENOME) == TRANSCRIPT_SEQ


@pytest.mark.parametrize("policy", ["off", GenomeMismatchPolicy.OFF])
def test_off(test_files, policy):
    fasta_filename, json_filename = test_files
    hdp = _hdp(json_filename, FastaSeqFetcher(fasta_filename, genome_mismatch=policy))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert hdp.seqfetcher.fetch_seq("NM_000314.4") == TRANSCRIPT_SEQ


def test_raise(test_files):
    fasta_filename, json_filename = test_files
    hdp = _hdp(json_filename, FastaSeqFetcher(fasta_filename, genome_mismatch="raise"))
    with pytest.raises(GenomeMismatchError, match="NM_000314.4"):
        hdp.seqfetcher.fetch_seq("NM_000314.4")
    assert issubclass(GenomeMismatchError, HGVSDataNotAvailableError)
    assert hdp.seqfetcher.fetch_seq(MATCHES_GENOME) == TRANSCRIPT_SEQ


def test_raise_falls_through_chained_seqfetcher(test_files):
    """ The point of raising: put FastaSeqFetcher first and fall back to a real transcript source """
    fasta_filename, json_filename = test_files
    seqfetcher = ChainedSeqFetcher(FastaSeqFetcher(fasta_filename, genome_mismatch="raise"), StaticSeqFetcher())
    hdp = _hdp(json_filename, seqfetcher)
    assert hdp.seqfetcher.fetch_seq("NM_000314.4") == "REAL"
    assert hdp.seqfetcher.fetch_seq(MATCHES_GENOME) == TRANSCRIPT_SEQ


def test_invalid_policy(test_files):
    fasta_filename, _ = test_files
    with pytest.raises(ValueError):
        FastaSeqFetcher(fasta_filename, genome_mismatch="loud")
