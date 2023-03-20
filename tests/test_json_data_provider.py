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
            print(tx_data)
            if tx_data["tx_ac"] == expected_transcript:
                found = True
                self.assertEqual(tx_data["alt_ac"], "NC_000007.13")
                continue
        self.assertTrue(found)

    def test_get_tx_for_region(self):
        found = False
        expected_transcript = "NM_001637.3"
        # Exonic coordinate
        for tx_data in self.json_data_provider.get_tx_for_region("NC_000007.13", "splign", 36570024, 36570025):
            if tx_data["tx_ac"] == expected_transcript:
                found = True
                self.assertEqual(tx_data["alt_strand"], -1)
                self.assertEqual(tx_data["start_i"], 36552548)
                self.assertEqual(tx_data["end_i"], 36764154)
                continue
        self.assertTrue(found)

    def test_get_tx_for_region_intron(self):
        """ Test case for https://github.com/SACGF/cdot/issues/38 """
        found = False
        expected_transcript = "NM_001637.3"
        # Coordinate below is intronic
        for tx_data in self.json_data_provider.get_tx_for_region("NC_000007.13", "splign", 36743533, 36745648):
            if tx_data["tx_ac"] == expected_transcript:
                found = True
                self.assertEqual(tx_data["alt_strand"], -1)
                self.assertEqual(tx_data["start_i"], 36552548)
                self.assertEqual(tx_data["end_i"], 36764154)
                continue
        self.assertTrue(found)


    def test_get_pro_ac_for_tx_ac(self):
        pro_ac = self.json_data_provider.get_pro_ac_for_tx_ac("NM_001637.3")
        self.assertEqual(pro_ac, "NP_001628.1")

    def test_get_gene_info(self):
        gene_info = self.json_data_provider.get_gene_info("GATA2")
        summary = gene_info.pop("summary")
        self.assertTrue("zinc-finger transcription factors" in summary)
        expected = {
            "hgnc": "GATA2",
            "maploc": "3q21.3",
            "descr": "GATA binding protein 2",
            "aliases": "{DCML,IMD21,MONOMAC,NFE1B}",
            "added": None,
        }
        self.assertEqual(gene_info, expected)


if __name__ == '__main__':
    unittest.main()
