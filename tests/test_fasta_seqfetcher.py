"""
Tests for FastaSeqFetcher warning/refusing to build a transcript that differs from the genome (#115)

Uses a tiny synthetic genome and transcripts. The genome_mismatch values are from the public examples
in #115: EGR2 NM_000399.3 (1 indel) and PTEN NM_000314.4 (1 substitution, 1 indel)
"""
import json
import warnings

from bioutils.sequences import reverse_complement
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


def test_fasta_opened_once(test_files, monkeypatch):
    """ FastaSeqFetcher shares one GenomeFastaSeqFetcher between contigs and exon-built transcripts (#169) """
    import cdot.hgvs.dataproviders.fasta_seqfetcher as fasta_seqfetcher_module

    fasta_filename, json_filename = test_files
    opened = []
    real_fasta_file = fasta_seqfetcher_module.FastaFile

    def counting_fasta_file(filename):
        opened.append(filename)
        return real_fasta_file(filename)

    monkeypatch.setattr(fasta_seqfetcher_module, "FastaFile", counting_fasta_file)
    seqfetcher = FastaSeqFetcher(fasta_filename)
    assert opened == [fasta_filename]
    assert seqfetcher.default_seqfetcher.contig_fastas is seqfetcher.contig_fastas

    hdp = _hdp(json_filename, seqfetcher)
    assert hdp.seqfetcher.fetch_seq(CONTIG, 10, 20) == GENOME[10:20]
    assert hdp.seqfetcher.fetch_seq(MATCHES_GENOME) == TRANSCRIPT_SEQ


def test_no_fasta_files():
    with pytest.raises(ValueError, match="at least one"):
        FastaSeqFetcher()
    with pytest.raises(ValueError, match="at least one"):
        ExonsFromGenomeFastaSeqFetcher()


# AAAAA run at 15-20, then 2 transcript bases not in the genome (filled with N)

PADDED_GENOME = "CGTCGTCGTC" + "CGTCG" + "AAAAA" + "CGTCGTCGTC" * 6
PADDED_TRANSCRIPTS = {
    "NM_000001.1": ("+", [[10, 30, 0, 1, 22, "M10 I2 M10"], [40, 60, 1, 23, 42, None]]),
    "NM_000002.1": ("-", [[10, 30, 0, 1, 22, "M5 I2 M15"]]),
    "NM_000003.1": ("+", [[10, 30, 0, 4, 23, None]]),
    "NM_000004.1": ("-", [[10, 30, 0, 4, 23, None]]),
    # Transcript bases 21-23 align nowhere, between the two exons
    "NM_000005.1": ("+", [[10, 30, 0, 1, 20, None], [40, 60, 1, 24, 43, None]]),
    "NM_000006.1": ("-", [[10, 30, 1, 24, 43, None], [40, 60, 0, 1, 20, None]]),
}


@pytest.fixture(scope="module")
def padded_files(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("fasta_padding")
    fasta_filename = str(tmp_path / "genome.fa")
    with open(fasta_filename, "w") as f:
        f.write(f">{CONTIG}\n{PADDED_GENOME}\n")
    pysam.faidx(fasta_filename)
    transcripts = {ac: {"id": ac, "gene_name": "TEST", "start_codon": 0, "stop_codon": exons[-1][4],
                        "genome_builds": {"GRCh37": {"contig": CONTIG, "strand": strand, "exons": exons}}}
                   for ac, (strand, exons) in PADDED_TRANSCRIPTS.items()}
    json_filename = str(tmp_path / "cdot.json")
    with open(json_filename, "w") as f:
        json.dump({"cdot_version": "0.2.36", "genome_builds": ["GRCh37"],
                   "transcripts": transcripts, "genes": {}}, f)
    return fasta_filename, json_filename


def _assembly_mapper(hdp):
    from hgvs.assemblymapper import AssemblyMapper
    return AssemblyMapper(hdp, assembly_name="GRCh37", alt_aln_method="splign", replace_reference=True)


@pytest.mark.parametrize("ac, positions", [
    ("NM_000001.1", "n.11_12"), ("NM_000002.1", "n.6_7"), ("NM_000003.1", "n.1_3"),
])
def test_slice_with_padding_raises(padded_files, ac, positions):
    fasta_filename, json_filename = padded_files
    hdp = _hdp(json_filename, FastaSeqFetcher(fasta_filename))
    with pytest.raises(GenomeMismatchError, match=rf"{ac} {positions} unknown"):
        hdp.seqfetcher.fetch_seq(ac, 0, 20)


def test_slice_without_padding(padded_files):
    fasta_filename, json_filename = padded_files
    hdp = _hdp(json_filename, FastaSeqFetcher(fasta_filename))
    assert hdp.seqfetcher.fetch_seq("NM_000001.1", 0, 10) == "CGTCGAAAAA"
    assert hdp.seqfetcher.fetch_seq("NM_000001.1", 12, 22) == "CGTCGTCGTC"


def test_whole_transcript_with_padding_warns(padded_files):
    fasta_filename, json_filename = padded_files
    hdp = _hdp(json_filename, FastaSeqFetcher(fasta_filename))
    with pytest.warns(GenomeMismatchWarning, match="NM_000001.1 n.11_12"):
        assert hdp.seqfetcher.fetch_seq("NM_000001.1") == "CGTCGAAAAANNCGTCGTCGTC" + PADDED_GENOME[40:60]


@pytest.mark.parametrize("genome_mismatch", ["off", "raise"])
def test_padding_ignores_genome_mismatch_policy(padded_files, genome_mismatch):
    fasta_filename, json_filename = padded_files
    hdp = _hdp(json_filename, FastaSeqFetcher(fasta_filename, genome_mismatch=genome_mismatch))
    with pytest.raises(GenomeMismatchError):
        hdp.seqfetcher.fetch_seq("NM_000001.1", 0, 20)
    with pytest.warns(GenomeMismatchWarning, match="NM_000001.1 n.11_12"):
        hdp.seqfetcher.fetch_seq("NM_000001.1")


def test_padding_falls_through_chained_seqfetcher(padded_files):
    fasta_filename, json_filename = padded_files
    hdp = _hdp(json_filename, ChainedSeqFetcher(FastaSeqFetcher(fasta_filename), StaticSeqFetcher()))
    assert hdp.seqfetcher.fetch_seq("NM_000001.1", 0, 20) == "REAL"


def test_indel_normalized_into_padding_raises(padded_files):
    """ Shifting the del 3' along the A run would stop at the Ns """
    import hgvs.parser
    fasta_filename, json_filename = padded_files
    parser = hgvs.parser.Parser()
    var_g = parser.parse("NC_000011.9:g.18del")

    am = _assembly_mapper(_hdp(json_filename, FastaSeqFetcher(fasta_filename)))
    with pytest.raises(GenomeMismatchError):
        am.g_to_c(var_g, "NM_000001.1")
    assert str(am.c_to_g(parser.parse("NM_000001.1:c.30G>A"))) == "NC_000011.9:g.48G>A"
    with pytest.raises(GenomeMismatchError):
        am.c_to_g(parser.parse("NM_000001.1:c.11A>G"))


@pytest.mark.parametrize("ac, n_pos, g_pos", [
    ("NM_000003.1", 4, 11), ("NM_000003.1", 23, 30),
    ("NM_000004.1", 4, 30), ("NM_000004.1", 23, 11),
    ("NM_000005.1", 20, 30), ("NM_000005.1", 24, 41), ("NM_000005.1", 43, 60),
    ("NM_000006.1", 20, 41), ("NM_000006.1", 24, 30), ("NM_000006.1", 43, 11),
])
def test_unaligned_transcript_bases_map(padded_files, ac, n_pos, g_pos):
    """ Unaligned transcript bases (n.1_3 at the start, n.21_23 between exons) are skipped """
    from hgvs.alignmentmapper import AlignmentMapper
    from hgvs.parser import Parser
    _, json_filename = padded_files
    hdp = _hdp(json_filename, StaticSeqFetcher())
    interval = Parser().parse_n_interval(str(n_pos))
    g_interval = AlignmentMapper(hdp, ac, CONTIG, "splign").n_to_g(interval)
    assert (g_interval.start.base, g_interval.end.base) == (g_pos, g_pos)


@pytest.mark.parametrize("ac, seq", [
    ("NM_000005.1", PADDED_GENOME[10:30] + "NNN" + PADDED_GENOME[40:60]),
    ("NM_000006.1", reverse_complement(PADDED_GENOME[40:60]) + "NNN" + reverse_complement(PADDED_GENOME[10:30])),
])
def test_unaligned_bases_between_exons_padded(padded_files, ac, seq):
    fasta_filename, json_filename = padded_files
    hdp = _hdp(json_filename, FastaSeqFetcher(fasta_filename))
    with pytest.warns(GenomeMismatchWarning, match=f"{ac} n.21_23"):
        assert hdp.seqfetcher.fetch_seq(ac) == seq
