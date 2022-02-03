import os
import unittest
from inspect import getsourcefile
from os.path import abspath

import pyhgvs

from cdot.pyhgvs.pyhgvs_transcript import JSONPyHGVSTranscriptFactory
from .genome import MockGenomeTestFile


class TestPyHGVS(unittest.TestCase):
    def test_transcript(self):
        this_file_dir = os.path.dirname(abspath(getsourcefile(lambda: 0)))
        test_json_file = os.path.join(this_file_dir, "test_data/cdot.refseq.grch37.json")
        factory = JSONPyHGVSTranscriptFactory([test_json_file])

        HGVS_C_TO_G = [
            ('NM_001637.3:c.1582G>A', 'NC_000007.13:g.36561662C>T'),
        ]

        genome = MockGenomeTestFile(
            db_filename='grch37.fa',
            filename=os.path.join(this_file_dir, 'test_data/grch37.genome'),
            create_data=False)

        for hgvs_c, expected_hgvs_g in HGVS_C_TO_G:
            result = pyhgvs.parse_hgvs_name(hgvs_c, genome, get_transcript=factory.get_transcript_grch37)
            name = pyhgvs.HGVSName(expected_hgvs_g)
            expected = (name.chrom, name.start, name.ref_allele, name.alt_allele)
            self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
