import os
import unittest
from inspect import getsourcefile
from os.path import abspath

from cdot.hgvs.dataproviders.json_data_provider import JSONDataProvider


class BatchRecordingJSONDataProvider(JSONDataProvider):
    """ Records how the transcript-walking methods retrieve their transcripts """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.batches = []
        self.single_fetches = []

    def _get_transcripts(self, tx_acs):
        tx_acs = list(tx_acs)
        self.batches.append(tx_acs)
        return {tx_ac: self._get_transcript(tx_ac) for tx_ac in tx_acs}

    def _get_transcript(self, tx_ac):
        self.single_fetches.append(tx_ac)
        return super()._get_transcript(tx_ac)


class TestRetrieveMany(unittest.TestCase):
    """ Every method that walks a list of transcripts must go through _get_transcripts, so that a
        provider able to answer for many accessions at once (DB query, bulk REST call) only has to
        override that one method - @see AbstractJSONDataProvider._get_transcripts """

    GENE = "AOAH"
    CONTIG = "NC_000007.13"

    def setUp(self):
        this_file_dir = os.path.dirname(abspath(getsourcefile(lambda: 0)))
        test_json_file = os.path.join(this_file_dir, "test_data/cdot.refseq.grch37.json")
        self.data_provider = BatchRecordingJSONDataProvider([test_json_file])

    def _assert_retrieved_in_batches(self, expected_batches=1):
        """ The point is that the number of retrievals is constant, not per-transcript """
        batches = self.data_provider.batches
        self.assertEqual(expected_batches, len(batches))
        for batch in batches:
            self.assertTrue(batch, "batch was not empty")

    def test_get_tx_for_gene_retrieves_in_one_batch(self):
        self.data_provider.get_tx_for_gene(self.GENE)
        self._assert_retrieved_in_batches()

    def test_get_tx_ac_tags_for_gene_retrieves_in_batches(self):
        # Two batches (rank by length, then tags), not one query per transcript - which
        # is the whole point: a few hundred transcripts cost 2 retrievals, not 400.
        self.data_provider.get_tx_ac_tags_for_gene(self.GENE, "GRCh37")
        self._assert_retrieved_in_batches(2)

    def test_get_tx_for_region_retrieves_in_one_batch(self):
        self.data_provider.get_tx_for_region(self.CONTIG, "splign", 36570024, 36570025)
        self._assert_retrieved_in_batches()

    def test_default_retrieves_one_at_a_time(self):
        """ Providers holding transcripts in memory keep working without an override """
        tx_acs = [tx_ac for tx_ac, _tags in self.data_provider.get_tx_ac_tags_for_gene(self.GENE, "GRCh37")]
        self.data_provider.single_fetches.clear()

        transcripts = JSONDataProvider._get_transcripts(self.data_provider, tx_acs)
        self.assertEqual(tx_acs, list(transcripts))
        self.assertEqual(tx_acs, self.data_provider.single_fetches)

    def test_duplicate_accessions_collapse(self):
        tx_acs = [tx_ac for tx_ac, _tags in self.data_provider.get_tx_ac_tags_for_gene(self.GENE, "GRCh37")]
        self.data_provider.single_fetches.clear()

        transcripts = JSONDataProvider._get_transcripts(self.data_provider, tx_acs + tx_acs)
        self.assertEqual(tx_acs, list(transcripts))
