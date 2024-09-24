import abc

from pysam.libcfaidx import FastaFile
from hgvs.dataproviders.interface import Interface
from hgvs.exceptions import HGVSDataNotAvailableError


class PrefixSeqFetcher:
    def __init__(self, default_seq_fetcher=None):
        self.default_seq_fetcher = default_seq_fetcher
        self.prefix_seq_fetchers = {}

    def add_seq_fetcher(self, prefix, seq_fetcher):
        self.prefix_seq_fetchers[prefix] = seq_fetcher

    def fetch_seq(self, ac, start_i=None, end_i=None):
        for prefix, sf in self.prefix_seq_fetchers.items():
            if ac.startswith(prefix):
                return sf.fetch_seq(ac, start_i=start_i, end_i=end_i)
        if self.default_seq_fetcher:
            return self.default_seq_fetcher.fetch_seq(ac, start_i=start_i, end_i=end_i)

        known_prefixes = ','.join(self.prefix_seq_fetchers.keys())
        msg = f"Couldn't handle '{ac}', must match known prefixes: '{known_prefixes}'. No default set"
        raise HGVSDataNotAvailableError(msg)


class ChainedSeqFetcher:
    """ This takes multiple SeqFetcher instances, and tries them in order if HGVSDataNotAvailableError
        until one succeeds (or finally throws)

        This is useful if you want to use FastaSeqFetcher (below) as a fallback if SeqFetcher fails

        seqfetcher = ChainedSeqFetcher(SeqFetcher(), FastaSeqFetcher(fasta_filename))
    """

    def __init__(self, *args):
        self.seq_fetchers = list(args)
        self.source = ", ".join(s.source for s in self.seq_fetchers)

    def set_data_provider(self, hdp: Interface):
        for seqfetcher in self.seq_fetchers:
            try:
                seqfetcher.set_data_provider(hdp)
            except AttributeError:
                pass

    def fetch_seq(self, ac, start_i=None, end_i=None):
        exceptions = []
        for sf in self.seq_fetchers:
            try:
                return sf.fetch_seq(ac, start_i=start_i, end_i=end_i)
            except HGVSDataNotAvailableError as e:
                exceptions.append(e)

        raise HGVSDataNotAvailableError(exceptions)


class AbstractTranscriptSeqFetcher:
    def __init__(self, *args, cache=True):
        self.cache = cache
        self.transcript_cache = {}
        self.hdp = None  # Set when passed to data provider (via set_data_provider)

    @abc.abstractmethod
    def _get_transcript_seq(self, ac, alt_ac, alt_aln_method):
        pass

    def get_transcript_seq(self, ac):
        transcript_seq = self.transcript_cache.get(ac)
        if not transcript_seq:
            transcript_seq = self._get_transcript_seq(ac)
            if self.cache:
                self.transcript_cache[ac] = transcript_seq
        return transcript_seq

    def set_data_provider(self, hdp: Interface):
        self.hdp = hdp

    def fetch_seq(self, ac, start_i=None, end_i=None):
        if self.hdp is None:
            raise HGVSDataNotAvailableError("You need to set set_data_provider() before calling fetch_seq()")

        transcript_seq = self.get_transcript_seq(ac)
        if start_i is None:
            start_i = 0
        if end_i is None:
            end_i = len(transcript_seq)
        return transcript_seq[start_i:end_i]
