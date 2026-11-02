""" cdot_json.add_gencode_hgnc - reads GENCODE transcript -> HGNC metadata into Ensembl data """
import gzip
import os
import shutil
import tempfile
import unittest

from generate_transcript_data.cdot_json import add_gencode_hgnc


class TestAddGencodeHGNC(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def _write(self, content: str) -> str:
        filename = os.path.join(self.tmp_dir, "gencode.metadata.HGNC.gz")
        with gzip.open(filename, "wt") as f:
            f.write(content)
        return filename

    @staticmethod
    def _genes_and_transcripts():
        genes = {"ENSG00000012048.23": {"gene_symbol": "BRCA1"}}
        transcripts = {"ENST00000357654.9": {"gene_symbol": "BRCA1", "gene_version": "ENSG00000012048.23"}}
        return genes, transcripts

    def test_blank_lines_skipped(self):
        """ #162 - a blank line used to raise ValueError """
        filename = self._write("\nENST00000357654.9\tBRCA1\tHGNC:1100\n\n")
        genes, transcripts = self._genes_and_transcripts()
        add_gencode_hgnc(filename, genes, transcripts)
        self.assertEqual("1100", transcripts["ENST00000357654.9"]["hgnc"])
        self.assertEqual("1100", genes["ENSG00000012048.23"]["hgnc"])

    def test_malformed_line(self):
        filename = self._write("ENST00000357654.9\tBRCA1\tHGNC:1100\textra\n")
        genes, transcripts = self._genes_and_transcripts()
        with self.assertRaisesRegex(ValueError, "Expected 2 or 3 TSV values"):
            add_gencode_hgnc(filename, genes, transcripts)
