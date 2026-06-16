"""
Core loading stores typed msgspec structs (cdot/models.py) instead of plain
dicts, for faster/smaller loading. The structs are dict-read-compatible so the
data-provider base methods - and external subclasses that return plain dicts
(e.g. the REST provider here, or cdot_rest's Redis provider) - keep working
unchanged.

These tests assert that a struct-backed store and a dict-backed store produce
*identical* public output, i.e. the typed storage is invisible to consumers.
"""
import json
import os
import unittest
from inspect import getsourcefile
from os.path import abspath

from cdot.hgvs.dataproviders.json_data_provider import JSONDataProvider

THIS_DIR = os.path.dirname(abspath(getsourcefile(lambda: 0)))
REFSEQ_JSON = os.path.join(THIS_DIR, "test_data/cdot.refseq.grch37.json")
ENSEMBL_JSON = os.path.join(THIS_DIR, "test_data/cdot.ensembl.grch38.json")
ALN = "splign"


class TestStructDictCompat(unittest.TestCase):
    def _assert_identical(self, json_file):
        # Struct-backed provider (the new default core loading)
        struct_dp = JSONDataProvider([json_file])

        # Dict-backed provider: same data, but the store holds plain dicts -
        # this is what a REST/Redis subclass returns from _get_transcript.
        dict_dp = JSONDataProvider([json_file])
        with open(json_file) as f:
            dict_dp.transcripts = json.load(f)["transcripts"]

        self.assertEqual(set(struct_dp.transcripts), set(dict_dp.transcripts))
        # confirm we really are comparing struct vs dict storage
        a_tx = next(iter(struct_dp.transcripts.values()))
        self.assertNotIsInstance(a_tx, dict)
        self.assertIsInstance(next(iter(dict_dp.transcripts.values())), dict)

        checked = 0
        for ac, tx in struct_dp.transcripts.items():
            for build in tx["genome_builds"].values():
                contig = build["contig"]
                self.assertEqual(
                    struct_dp.get_tx_exons(ac, contig, ALN),
                    dict_dp.get_tx_exons(ac, contig, ALN),
                )
                self.assertEqual(
                    struct_dp.get_tx_info(ac, contig, ALN),
                    dict_dp.get_tx_info(ac, contig, ALN),
                )
                self.assertEqual(
                    struct_dp.get_tx_mapping_options(ac),
                    dict_dp.get_tx_mapping_options(ac),
                )
                checked += 1
            self.assertEqual(
                struct_dp.get_tx_identity_info(ac),
                dict_dp.get_tx_identity_info(ac),
            )
            self.assertEqual(
                struct_dp.get_pro_ac_for_tx_ac(ac),
                dict_dp.get_pro_ac_for_tx_ac(ac),
            )
        self.assertGreater(checked, 0)

    def test_refseq_struct_matches_dict(self):
        self._assert_identical(REFSEQ_JSON)

    def test_ensembl_struct_matches_dict(self):
        self._assert_identical(ENSEMBL_JSON)

    def test_get_tx_for_region_matches(self):
        struct_dp = JSONDataProvider([REFSEQ_JSON])
        dict_dp = JSONDataProvider([REFSEQ_JSON])
        with open(REFSEQ_JSON) as f:
            dict_dp.transcripts = json.load(f)["transcripts"]

        # region spanning the AOAH transcript on GRCh37
        contig = "NC_000007.13"
        self.assertEqual(
            struct_dp.get_tx_for_region(contig, ALN, 36552548, 36764154),
            dict_dp.get_tx_for_region(contig, ALN, 36552548, 36764154),
        )


if __name__ == "__main__":
    unittest.main()
