import os
import unittest
from inspect import getsourcefile
from os.path import abspath

import pyhgvs

from cdot.pyhgvs.pyhgvs_transcript import JSONPyHGVSTranscriptFactory, is_sacgf_pyhgvs_fork
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

        sacgf_pyhgvs_fork = is_sacgf_pyhgvs_fork()

        def get_transcript(transcript_id):
            return factory.get_transcript_grch37(transcript_id, sacgf_pyhgvs_fork=sacgf_pyhgvs_fork)

        for hgvs_c, expected_hgvs_g in HGVS_C_TO_G:
            result = pyhgvs.parse_hgvs_name(hgvs_c, genome, get_transcript=get_transcript)
            name = pyhgvs.HGVSName(expected_hgvs_g)
            expected = (name.chrom, name.start, name.ref_allele, name.alt_allele)
            self.assertEqual(result, expected)

    def test_non_coding_transcript(self):
        this_file_dir = os.path.dirname(abspath(getsourcefile(lambda: 0)))
        test_json_file = os.path.join(this_file_dir, "test_data/cdot.refseq.grch37.json")
        factory = JSONPyHGVSTranscriptFactory([test_json_file])

        genome = MockGenomeTestFile(
            db_filename='grch37.fa',
            filename=os.path.join(this_file_dir, 'test_data/grch37.genome'),
            create_data=False)

        transcript_id = "NR_023343.1"
        sacgf_pyhgvs_fork = is_sacgf_pyhgvs_fork()
        pyhgvs_transcript = factory.get_transcript_grch37(transcript_id, sacgf_pyhgvs_fork=sacgf_pyhgvs_fork)
        self.assertFalse(pyhgvs_transcript.is_coding, f"Transcript {transcript_id} is non-coding")


if __name__ == '__main__':
    unittest.main()
