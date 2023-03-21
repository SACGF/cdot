import re

from pysam.libcfaidx import FastaFile
from hgvs.dataproviders.interface import Interface
from hgvs.exceptions import HGVSDataNotAvailableError
from bioutils.sequences import reverse_complement


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


class FastaSeqFetcher:
    """ This produces artificial transcript sequences by pasting together exons from the genome
        It is possible that this does not exactly match the transcript sequences - USE AT OWN RISK! """
    def __init__(self, *args, cache=True):
        self.cache = cache
        self.transcript_cache = {}
        self.hdp = None  # Set when passed to data provider (via set_data_provider)
        self.source = "Local Fasta file reference"
        self.contig_fastas = {}
        self.cigar_pattern = re.compile(r"(\d+)([=DIX])")
        for fasta_filename in args:
            fasta_file = FastaFile(fasta_filename)
            for contig in fasta_file.references:
                self.contig_fastas[contig] = fasta_file

        if not self.contig_fastas:
            raise ValueError("Need to provide at least one of fasta file as argument")

    def set_data_provider(self, hdp: Interface):
        self.hdp = hdp

    def fetch_seq(self, ac, start_i=None, end_i=None):
        if fasta_file := self.contig_fastas.get(ac):  # Contig
            return fasta_file.fetch(ac, start_i, end_i).upper()

        if self.hdp is None:
            raise HGVSDataNotAvailableError("You need to set set_data_provider() before calling fetch_seq()")

        possible_contigs = set()
        for tx_mo in self.hdp.get_tx_mapping_options(ac):
            alt_ac = tx_mo["alt_ac"]
            possible_contigs.add(alt_ac)
            if alt_ac in self.contig_fastas:
                transcript_seq = self._get_transcript_seq(ac, alt_ac, tx_mo["alt_aln_method"])
                if start_i is None:
                    start_i = 0
                if end_i is None:
                    end_i = len(transcript_seq)
                return transcript_seq[start_i:end_i]

        msg = f"Failed to fetch {ac} from {self.source}. "
        if possible_contigs:
            possible_contigs = sorted(possible_contigs)
            raise HGVSDataNotAvailableError(f"{msg} No Fasta provided with contigs: {possible_contigs}")
        raise HGVSDataNotAvailableError(f"{msg} Transcript '{ac}' not found.")

    def _get_transcript_seq(self, ac, alt_ac, alt_aln_method):
        transcript_seq = self.transcript_cache.get(ac)
        if not transcript_seq:
            transcript_seq = self._fetch_seq_from_fasta(ac, alt_ac, alt_aln_method)
            if self.cache:
                self.transcript_cache[ac] = transcript_seq
        return transcript_seq

    def _fetch_seq_from_fasta(self, ac, alt_ac, alt_aln_method):
        fasta_file = self.contig_fastas[alt_ac]

        exons = self.hdp.get_tx_exons(ac, alt_ac, alt_aln_method)
        exon_sequences = []
        for exon in sorted(exons, key=lambda ex: ex["ord"]):
            exon_seq = fasta_file.fetch(alt_ac, exon["alt_start_i"], exon["alt_end_i"])
            exon_seq = exon_seq.upper()

            # Cigar is mapping from transcript to genome.
            # We are going from genome to transcript so operations are reversed
            exon_seq_list = []
            start = 0
            for (length_str, op) in self.cigar_pattern.findall(exon["cigar"]):
                length = int(length_str)
                if op == 'D':
                    exon_seq_list.append("N" * length)
                elif op == 'I':
                    pass  # leave out
                else:  # match/mismatch
                    exon_seq_list.append(exon_seq[start:start+length])
                    start += length

            exon_seq = "".join(exon_seq_list)
            if exon["alt_strand"] == -1:
                exon_seq = reverse_complement(exon_seq)
            exon_sequences.append(exon_seq)
        return "".join(exon_sequences)
