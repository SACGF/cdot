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
        test_json_file = os.path.join(this_file_dir, "test_data/cdot.refseq.grch37.json")
        test_transcripts_file = os.path.join(this_file_dir, "test_data/transcript_sequences.json")
        mock_seqfetcher = MockSeqFetcher(test_transcripts_file)
        seqfetcher = ChainedSeqFetcher(mock_seqfetcher, SeqFetcher())
        cls.json_data_provider = JSONDataProvider([test_json_file], seqfetcher=seqfetcher)

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

    def test_get_tx_info(self):
        # We only have data for GRCh37 but none for 38

        # Make sure 37 works
        tx_info = self.json_data_provider.get_tx_info("NM_001637.3", "NC_000007.13", "splign")
        print(tx_info)

        # Make sure 38 fails
        with self.assertRaises(HGVSDataNotAvailableError):
            tx_info = self.json_data_provider.get_tx_info("NM_001637.3", "NC_000007.14", "splign")

    def test_get_tx_exons_unknown_contig(self):
        """ A contig outside the provider's assemblies is missing data, the error biocommons
            hgvs handles (the same one get_tx_info raises), not a ValueError """
        with self.assertRaisesRegex(HGVSDataNotAvailableError, "Contig 'chr7' not supported"):
            self.json_data_provider.get_tx_exons("NM_001637.3", "chr7", "splign")


class VersionlessIdJSONDataProvider(JSONDataProvider):
    """ Keys transcripts by versionless accession, like a store that holds one version per
        accession (eg cdot_rest's Redis provider), so the id in hand is "NM_001637" while the
        record's "id" is "NM_001637.3" """

    @staticmethod
    def _versionless(tx_ac):
        return tx_ac.rpartition(".")[0] or tx_ac

    def _get_transcript(self, tx_ac):
        if transcript := super()._get_transcript(tx_ac):
            return transcript
        for versioned_tx_ac, transcript in self.transcripts.items():
            if self._versionless(versioned_tx_ac) == tx_ac:
                return transcript
        return None

    def _get_transcript_ids_for_gene(self, gene):
        return {self._versionless(tx_ac) for tx_ac in super()._get_transcript_ids_for_gene(gene)}

    def _get_contig_interval_tree(self, alt_ac):
        from intervaltree import IntervalTree
        tree = super()._get_contig_interval_tree(alt_ac)
        return IntervalTree(iv._replace(data=self._versionless(iv.data)) for iv in tree)


class TestVersionlessIdProvider(unittest.TestCase):
    """ Every transcript-listing method returns the record's versioned accession as tx_ac, whatever
        id the provider stores its transcripts under """

    @classmethod
    def setUpClass(cls):
        this_file_dir = os.path.dirname(abspath(getsourcefile(lambda: 0)))
        test_json_file = os.path.join(this_file_dir, "test_data/cdot.refseq.grch37.json")
        cls.json_data_provider = VersionlessIdJSONDataProvider([test_json_file])

    def test_get_tx_for_region_returns_versioned_accession(self):
        tx_acs = {tx["tx_ac"] for tx in
                  self.json_data_provider.get_tx_for_region("NC_000007.13", "splign", 36570024, 36570025)}
        self.assertIn("NM_001637.3", tx_acs)
        self.assertNotIn("NM_001637", tx_acs)

    def test_get_tx_for_gene_returns_versioned_accession(self):
        tx_acs = {tx["tx_ac"] for tx in self.json_data_provider.get_tx_for_gene("AOAH")}
        self.assertIn("NM_001637.3", tx_acs)
        self.assertNotIn("NM_001637", tx_acs)

    def test_get_tx_ac_tags_for_gene_returns_versioned_accession(self):
        tx_acs = {tx_ac for tx_ac, _tags in self.json_data_provider.get_tx_ac_tags_for_gene("AOAH", "GRCh37")}
        self.assertIn("NM_001637.3", tx_acs)
        self.assertNotIn("NM_001637", tx_acs)


if __name__ == '__main__':
    unittest.main()
