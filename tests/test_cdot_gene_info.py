""" cdot_gene_info.py (Entrez calls are mocked, these tests do not use the network) """
import gzip
import json
import os
import shutil
import sys
import tempfile
import unittest
from datetime import datetime
from unittest.mock import patch

from generate_transcript_data import cdot_gene_info
from generate_transcript_data.json_schema_version import JSON_SCHEMA_VERSION


class FakeGeneSummary(dict):
    """ Entrez.read returns dicts with an extra 'attributes' """
    def __init__(self, uid, **kwargs):
        super().__init__(**kwargs)
        self.attributes = {"uid": uid}


class TestCdotGeneInfo(unittest.TestCase):
    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp_dir)

    def test_entrez_retries_exhausted(self):
        """ #161 - used to return None which the caller then tried to iterate """
        with patch.object(cdot_gene_info.Entrez, "epost", side_effect=IOError("NCBI down")) as mock_epost, \
                patch.object(cdot_gene_info.logging, "warning"):
            with self.assertRaisesRegex(RuntimeError, r"3 attempts for batch of 2 gene IDs \(672\.\.675\)") as cm:
                cdot_gene_info._get_entrez_gene_summary(["672", "675"])
        self.assertEqual(3, mock_epost.call_count)
        self.assertIsInstance(cm.exception.__cause__, IOError)

    def test_main_output(self):
        gene_info_filename = os.path.join(self.tmp_dir, "Homo_sapiens.gene_info.gz")
        with gzip.open(gene_info_filename, "wt") as f:
            f.write("#tax_id\tGeneID\tSymbol_from_nomenclature_authority\n")
            f.write("9606\t672\tBRCA1\n")
            f.write("9606\t999999\t-\n")  # No symbol, skipped
        mtime = datetime(2024, 1, 2, 3, 4, 5).timestamp()
        os.utime(gene_info_filename, (mtime, mtime))
        output_filename = os.path.join(self.tmp_dir, "gene_info.json.gz")

        summaries = [FakeGeneSummary("672", NomenclatureSymbol="BRCA1", MapLocation="17q21.31",
                                     NomenclatureName="BRCA1 DNA repair associated",
                                     OtherAliases="BRCAI", Summary="Tumour suppressor")]
        argv = ["cdot_gene_info.py", "--email", "test@example.com",
                "--gene-info", gene_info_filename, "--output", output_filename]
        with patch.object(sys, "argv", argv), \
                patch.object(cdot_gene_info, "_get_entrez_gene_summary", return_value=summaries) as mock_summary, \
                patch("builtins.print"):
            cdot_gene_info.main()
        mock_summary.assert_called_once_with(["672"])

        with gzip.open(output_filename, "rt") as f:
            data = json.load(f)
        self.assertEqual(JSON_SCHEMA_VERSION, data["cdot_version"])
        self.assertEqual("2024-01-02T03:04:05", data["gene_info_date"])
        self.assertEqual("BRCA1", data["gene_info"]["672"]["gene_symbol"])
