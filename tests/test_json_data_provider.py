import os
import unittest
from inspect import getsourcefile
from os.path import abspath

import hgvs
from hgvs.assemblymapper import AssemblyMapper

from cdot.hgvs.dataproviders.json_data_provider import JSONDataProvider


class TestJSONDataProvider(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        this_file_dir = os.path.dirname(abspath(getsourcefile(lambda: 0)))
        #        parent_dir = os.path.dirname(this_file_dir)
        test_json_file = os.path.join(this_file_dir, "test_data/cdot.refseq.grch37.json")
        cls.json_data_provider = JSONDataProvider([test_json_file])

    def test_transcript(self):
        am = AssemblyMapper(self.json_data_provider,
                            assembly_name='GRCh37', alt_aln_method='splign', replace_reference=True)
        HGVS_C_TO_G = [
            ('NM_001637.3:c.1582G>A', 'NC_000007.13:g.36561662C>T'),
        ]

        hp = hgvs.parser.Parser()
        for hgvs_c, expected_hgvs_g in HGVS_C_TO_G:
            var_c = hp.parse_hgvs_variant(hgvs_c)
            var_g = am.c_to_g(var_c)
            self.assertEqual(str(var_g), expected_hgvs_g)

    def test_get_tx_for_gene(self):
        found = False
        expected_transcript = "NM_001637.3"
        for tx_data in self.json_data_provider.get_tx_for_gene("AOAH"):
            (_gene, _cds_start_i, _cds_end_i, transcript_id, contig, _method) = tx_data
            if transcript_id == expected_transcript:
                found = True
                self.assertEqual(contig, "NC_000007.13")
                continue
        self.assertTrue(found)

    def test_get_tx_for_region(self):
        found = False
        expected_transcript = "NM_001637.3"
        for tx_data in self.json_data_provider.get_tx_for_region("NC_000007.13", "splign", 36570024, 36570025):
            (transcript_id, _alt_ac, strand, _method, tx_start, tx_end) = tx_data
            if transcript_id == expected_transcript:
                found = True
                self.assertEqual(strand, -1)
                self.assertEqual(tx_start, 36552548)
                self.assertEqual(tx_end, 36764154)
                continue
        self.assertTrue(found)



if __name__ == '__main__':
    unittest.main()
