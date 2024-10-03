import abc

from pysam.libcfaidx import FastaFile
from hgvs.dataproviders.interface import Interface
from hgvs.exceptions import HGVSDataNotAvailableError


class PrefixSeqFetcher:
    def __init__(self, default_seqfetcher=None):
        self.default_seqfetcher = default_seqfetcher
        self.prefix_seqfetchers = {}

    def add_seqfetcher(self, prefix, seqfetcher):
        self.prefix_seqfetchers[prefix] = seqfetcher

    @property
    def all_seqfetchers(self):
        seqfetchers = list(self.prefix_seqfetchers.values())
        if self.default_seqfetcher:
            seqfetchers.append(self.default_seqfetcher)
        return seqfetchers

    def set_data_provider(self, hdp: Interface):
        for seqfetcher in self.all_seqfetchers:
            try:
                seqfetcher.set_data_provider(hdp)
            except AttributeError:
                pass

    def fetch_seq(self, ac, start_i=None, end_i=None):
        for prefix, sf in self.prefix_seqfetchers.items():
            if ac.startswith(prefix):
                return sf.fetch_seq(ac, start_i=start_i, end_i=end_i)
        if self.default_seqfetcher:
            return self.default_seqfetcher.fetch_seq(ac, start_i=start_i, end_i=end_i)

        known_prefixes = ','.join(self.prefix_seqfetchers.keys())
        msg = f"Couldn't handle '{ac}', must match known prefixes: '{known_prefixes}'. No default set"
        raise HGVSDataNotAvailableError(msg)


class ChainedSeqFetcher:
    """ This takes multiple SeqFetcher instances, and tries them in order if HGVSDataNotAvailableError
        until one succeeds (or finally throws)

        This is useful if you want to use FastaSeqFetcher (below) as a fallback if SeqFetcher fails

        seqfetcher = ChainedSeqFetcher(SeqFetcher(), FastaSeqFetcher(fasta_filename))
    """

    def __init__(self, *args):
        self.seqfetchers = list(args)

    def set_data_provider(self, hdp: Interface):
        for seqfetcher in self.seqfetchers:
            try:
                seqfetcher.set_data_provider(hdp)
            except AttributeError:
                pass

    def fetch_seq(self, ac, start_i=None, end_i=None):
        exceptions = []
        for sf in self.seqfetchers:
            try:
                return sf.fetch_seq(ac, start_i=start_i, end_i=end_i)
            except HGVSDataNotAvailableError as e:
                exceptions.append(e)

        raise HGVSDataNotAvailableError(exceptions)

    @property
    def source(self):
        # This needs to execute after set_data_provider is called
        return ", ".join(s.source for s in self.seqfetchers)



class AbstractTranscriptSeqFetcher:
    def __init__(self, *args, cache=True):
        self.cache = cache
        self.transcript_cache = {}
        self.hdp = None  # Set when passed to data provider (via set_data_provider)

    @abc.abstractmethod
    def _get_transcript_seq(self, ac):
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
            raise HGVSDataNotAvailableError(f"{self}: You need to set set_data_provider() before calling fetch_seq()")

        transcript_seq = self.get_transcript_seq(ac)
        if start_i is None:
            start_i = 0
        if end_i is None:
            end_i = len(transcript_seq)
        return transcript_seq[start_i:end_i]
