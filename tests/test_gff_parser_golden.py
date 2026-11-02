""" Golden-file tests for the GTF/GFF3 parsers (issue #101)

    Every GTF/GFF3 fixture in tests/test_data is parsed and the genes/transcripts compared, byte for
    byte, against the JSON checked in under tests/test_data/gff_parser_golden/. The parser refactor
    must not change any of these.

    To regenerate after an intentional change to the output:

        CDOT_UPDATE_GOLDEN=1 python -m pytest tests/test_gff_parser_golden.py

    then review the diff of tests/test_data/gff_parser_golden/ before committing.
"""
import json
import logging
import os
import unittest
from inspect import getsourcefile

from generate_transcript_data.gff_parser import GTFParser, GFF3Parser
from generate_transcript_data.json_encoders import SortedSetEncoder

THIS_FILE_DIR = os.path.dirname(os.path.abspath(getsourcefile(lambda: 0)))
TEST_DATA_DIR = os.path.join(THIS_FILE_DIR, "test_data")
GOLDEN_DIR = os.path.join(TEST_DATA_DIR, "gff_parser_golden")
FAKE_URL = "http://fake.url"

# (fixture filename, parser class, genome build, parser kwargs[, golden suffix])
# The golden suffix lets the same fixture be parsed again with different options and compared to its own golden
FIXTURES = [
    # Ensembl GTF
    ("ensembl_test.GRCh38.104.gtf", GTFParser, "GRCh38", {}),
    ("ensembl_test.GRCh38.111.gtf", GTFParser, "GRCh38", {}),
    ("ensembl_test.GRCh38.108.selenoprotein.gtf", GTFParser, "GRCh38", {}),
    ("ensembl_test.GRCh38.115.selenoprotein.gtf", GTFParser, "GRCh38", {}),
    ("ensembl_test.GRCh38.115.MT.gtf", GTFParser, "GRCh38", {}),
    ("ensembl_test.GRCh38.116.chr21_slice.gtf.gz", GTFParser, "GRCh38", {}),
    # Ensembl GFF3 (protein versions on CDS rows from release 114)
    ("ensembl_test.GRCh38.115.MT.gff3", GFF3Parser, "GRCh38", {}),
    # UCSC GTF (no versions on anything)
    ("hg19_chrY_300kb_genes.gtf", GTFParser, "GRCh37", {}),
    # RefSeq GFF3
    ("refseq_test.GRCh38.p13_genomic.109.20210514.gff", GFF3Parser, "GRCh38", {}),
    ("refseq_test.GRCh38.p14_genomic.RS_2023_03.gff", GFF3Parser, "GRCh38", {}),
    ("refseq_test.RS_2025_08.chr21_slice.gff.gz", GFF3Parser, "GRCh38", {}),
    ("refseq_test.transcript_coordinate_hole.gff", GFF3Parser, "GRCh38", {}),
    ("refseq_test.historical_RS_2024_08.gff", GFF3Parser, "GRCh38", {"skip_missing_parents": True}),
    ("refseq_grch37_mt.gff", GFF3Parser, "GRCh37", {}),
    ("refseq_grch38.p14_mt.gff", GFF3Parser, "GRCh38", {}),
    ("refseq_test.selenoprotein.gff", GFF3Parser, "GRCh38", {}),
    ("refseq_test.RS_2024_08.selenoprotein.gff", GFF3Parser, "GRCh38", {}),
    ("refseq_test.RS_2025_08.ribosomal_slippage.gff", GFF3Parser, "GRCh38", {}),
    ("refseq_test.GRCh37.105.20220307.ACTN3.gff", GFF3Parser, "GRCh37", {}),
    ("refseq_test.RS_2025_08.ACTN3.gff", GFF3Parser, "GRCh38", {}),
    # Synthetic: PAR gene on X and Y (other_chroms), CDS starting outside the exons (codons_unplaced),
    # gene on an NT_ scaffold (discarded by default). Parsed again keeping that scaffold and without
    # contig conversion, for a genome build that has no bioutils assembly
    ("refseq_test.par_and_unplaced_codon.gff", GFF3Parser, "GRCh38", {}),
    ("refseq_test.par_and_unplaced_codon.gff", GFF3Parser, "CustomBuild",
     {"no_contig_conversion": True, "discard_contigs_with_underscores": False}, "no_contig_conversion"),
    # Synthetic: Ensembl GRCh37 GTFs repeat PAR genes on Y with the same IDs (other_chroms)
    ("ensembl_test.GRCh37.87.par.gtf", GTFParser, "GRCh37", {}),
]


def _golden_filename(fixture, golden_suffix=""):
    if golden_suffix:
        golden_suffix = "." + golden_suffix
    return os.path.join(GOLDEN_DIR, fixture.removesuffix(".gz") + golden_suffix + ".json")


def _parse(fixture, parser_class, genome_build, kwargs) -> str:
    parser = parser_class(os.path.join(TEST_DATA_DIR, fixture), genome_build, FAKE_URL, **kwargs)
    logging.disable(logging.WARNING)  # Some fixtures deliberately have unplaceable features
    try:
        genes, transcripts = parser.get_genes_and_transcripts()
    finally:
        logging.disable(logging.NOTSET)
    data = {"genes": genes, "transcripts": transcripts}
    return json.dumps(data, cls=SortedSetEncoder, sort_keys=True, indent=1) + "\n"


class TestGFFParserGolden(unittest.TestCase):
    def test_all_fixtures_have_golden(self):
        fixtures = {f for f in os.listdir(TEST_DATA_DIR)
                    if f.endswith((".gtf", ".gtf.gz", ".gff", ".gff.gz", ".gff3", ".gff3.gz"))}
        covered = {f[0] for f in FIXTURES}
        self.assertEqual(fixtures - covered, set(), "GTF/GFF fixtures without a golden file")


def _make_test(fixture, parser_class, genome_build, kwargs, golden_suffix=""):
    def test(self):
        actual = _parse(fixture, parser_class, genome_build, kwargs)
        golden_filename = _golden_filename(fixture, golden_suffix)
        if os.environ.get("CDOT_UPDATE_GOLDEN"):
            os.makedirs(GOLDEN_DIR, exist_ok=True)
            with open(golden_filename, "w") as f:
                f.write(actual)
        with open(golden_filename) as f:
            expected = f.read()
        self.assertEqual(expected, actual, f"{fixture} differs from {os.path.relpath(golden_filename)}")
    return test


for _fixture, _parser_class, _genome_build, _kwargs, *_golden_suffix in FIXTURES:
    _name = "test_" + "_".join([_fixture] + _golden_suffix).replace(".", "_")
    setattr(TestGFFParserGolden, _name, _make_test(_fixture, _parser_class, _genome_build, _kwargs, *_golden_suffix))
