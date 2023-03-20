from pysam.libcfaidx import FastaFile
from hgvs.dataproviders.interface import Interface
from bioutils.sequences import reverse_complement


class FastaSeqFetcher:
    """
        This produces artificial transcript sequences by pasting together exons from the genome
        It is possible that this does not exactly match the transcript sequences - USE AT OWN RISK!
    """
    def __init__(self, *args):
        self.hdp = None  # Set when passed to data provider (via set_data_provider)
        self.source = "Local Fasta file reference"
        self.contig_fastas = {}
        for fasta_filename in args:
            fasta_file = FastaFile(fasta_filename)
            for contig in fasta_file.references:
                self.contig_fastas[contig] = fasta_file

        if not self.contig_fastas:
            raise ValueError("Need to provide at least one of fasta file as argument")

    def set_data_provider(self, hdp: Interface):
        self.hdp = hdp

    def fetch_seq(self, ac, start_i=None, end_i=None):
        if self.hdp is None:
            raise HGVSDataNotAvailableError("You need to set set_data_provider() before calling fetch_seq()")

        possible_contigs = set()
        for tx_mo in self.hdp.get_tx_mapping_options(ac):
            alt_ac = tx_mo["alt_ac"]
            possible_contigs.add(alt_ac)
            if alt_ac in self.contig_fastas:
                return self._fetch_seq_from_fasta(ac, alt_ac, tx_mo["alt_aln_method"], start_i=start_i, end_i=end_i)

        possible_contigs = sorted(possible_contigs)
        raise HGVSDataNotAvailableError(f"Failed to fetch {ac} from {self.source}. "
                                        f"No Fasta provided with contigs: {possible_contigs}")

    def _fetch_seq_from_fasta(self, ac, alt_ac, alt_aln_method, start_i=None, end_i=None):
        # TODO: This doesn't yet handle gaps
        # TODO: Doesn't handle start_id or end_i being non-None

        fasta_file = self.contig_fastas[alt_ac]

        exons = self.hdp.get_tx_exons(ac, alt_ac, alt_aln_method)
        exon_sequences = []
        for exon in sorted(exons, key=lambda ex: ex["ord"]):
            seq = fasta_file.fetch(alt_ac, exon["alt_start_i"], exon["alt_end_i"])
            if exon["alt_strand"] == -1:
                seq = reverse_complement(seq)
            exon_sequences.append(seq)

        return "".join(exon_sequences)
