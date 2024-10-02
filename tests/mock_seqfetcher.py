import json

from hgvs.exceptions import HGVSDataNotAvailableError


class MockSeqFetcher:
    def __init__(self, filename):
        with open(filename) as f:
            self.transcripts = json.load(f)
        self.source = f"Mock: Local JSON file: {filename}"

    def fetch_seq(self, ac, start_i=None, end_i=None):
        seq = self.transcripts.get(ac)
        if seq is None:
            raise HGVSDataNotAvailableError()
        if start_i is None:
            start_i = 0
        if end_i is None:
            end_i = len(seq)
        return seq[start_i:end_i]

