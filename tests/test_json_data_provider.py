import os
import unittest
from inspect import getsourcefile
from os.path import abspath

import hgvs
from hgvs.assemblymapper import AssemblyMapper

from cdot.hgvs.dataproviders.json_data_provider import JSONDataProvider


class TestJSONDataProvider(unittest.TestCase):
    def test_transcript(self):
        this_file_dir = os.path.dirname(abspath(getsourcefile(lambda: 0)))
#        parent_dir = os.path.dirname(this_file_dir)
        test_json_file = os.path.join(this_file_dir, "test_data/cdot.refseq.grch37.json")
        json_data_provider = JSONDataProvider(assembly_json={"GRCh37": test_json_file})
        am = AssemblyMapper(json_data_provider,
                            assembly_name='GRCh37', alt_aln_method='splign', replace_reference=True)
        HGVS_C_TO_G = [
            ('NM_001637.3:c.1582G>A', 'NC_000007.13:g.36561662C>T'),
        ]

        hp = hgvs.parser.Parser()
        for hgvs_c, expected_hgvs_g in HGVS_C_TO_G:
            var_c = hp.parse_hgvs_variant(hgvs_c)
            var_g = am.c_to_g(var_c)
            self.assertEqual(str(var_g), expected_hgvs_g)


if __name__ == '__main__':
    unittest.main()
