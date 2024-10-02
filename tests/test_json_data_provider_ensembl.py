import os
import unittest
from inspect import getsourcefile
from os.path import abspath

import hgvs
from hgvs.assemblymapper import AssemblyMapper
from hgvs.dataproviders.seqfetcher import SeqFetcher
from hgvs.exceptions import HGVSDataNotAvailableError

from cdot.hgvs.dataproviders import ChainedSeqFetcher
from cdot.hgvs.dataproviders.json_data_provider import JSONDataProvider
from tests.mock_seqfetcher import MockSeqFetcher


class TestJSONDataProvider(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        this_file_dir = os.path.dirname(abspath(getsourcefile(lambda: 0)))
        #        parent_dir = os.path.dirname(this_file_dir)
        test_json_file = os.path.join(this_file_dir, "test_data/cdot.ensembl.grch38.json")
        test_transcripts_file = os.path.join(this_file_dir, "test_data/transcript_sequences.json")
        mock_seqfetcher = MockSeqFetcher(test_transcripts_file)
        seqfetcher = ChainedSeqFetcher(mock_seqfetcher, SeqFetcher())
        cls.json_data_provider = JSONDataProvider([test_json_file], seqfetcher=seqfetcher)

    def test_transcript(self):
        am = AssemblyMapper(self.json_data_provider,
                            assembly_name='GRCh38', alt_aln_method='splign', replace_reference=True)
        HGVS_C_TO_G = [
            ('ENST00000617537.5:c.1582G>A', 'NC_000007.14:g.36522056C>T'),
        ]

        hp = hgvs.parser.Parser()
        for hgvs_c, expected_hgvs_g in HGVS_C_TO_G:
            var_c = hp.parse_hgvs_variant(hgvs_c)
            var_g = am.c_to_g(var_c)
            self.assertEqual(str(var_g), expected_hgvs_g)

    def test_get_tx_for_gene(self):
        found = False
        expected_transcript = "ENST00000617537.5"
        for tx_data in self.json_data_provider.get_tx_for_gene("AOAH"):
            print(tx_data)
            if tx_data["tx_ac"] == expected_transcript:
                found = True
                self.assertEqual(tx_data["alt_ac"], "NC_000007.14")
                continue
        self.assertTrue(found)

    def test_get_tx_for_region(self):
        found = False
        expected_transcript = "ENST00000617537.5"
        # Exonic coordinate
        for tx_data in self.json_data_provider.get_tx_for_region("NC_000007.14", "splign", 36530417, 36530514):
            if tx_data["tx_ac"] == expected_transcript:
                found = True
                self.assertEqual(tx_data["alt_strand"], -1)
                self.assertEqual(tx_data["start_i"], 36512940)
                self.assertEqual(tx_data["end_i"], 36724494)
                continue
        self.assertTrue(found)

    def test_get_pro_ac_for_tx_ac(self):
        pro_ac = self.json_data_provider.get_pro_ac_for_tx_ac("ENST00000617537.5")
        self.assertEqual(pro_ac, "ENSP00000483783.1")

    def test_get_tx_info(self):
        # We only have data for GRCh38 but none for 37

        # Make sure 37 fails
        with self.assertRaises(HGVSDataNotAvailableError):
            tx_info = self.json_data_provider.get_tx_info("ENST00000617537.5", "NC_000007.13", "splign")

        # Make sure 38 works
        tx_info = self.json_data_provider.get_tx_info("ENST00000617537.5", "NC_000007.14", "splign")
        print(tx_info)


if __name__ == '__main__':
    unittest.main()
